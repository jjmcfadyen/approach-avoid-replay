library(lme4)
library(lmerTest)
library(interactions)
library(car)
library(emmeans)
library(tidyverse)
library(viridis)
library(extrafont)
library(gridExtra)
library(MuMIn)

#######################################
# LOAD DATA
#######################################

d <- read.csv("D:/2020_RiskyReplay/results/replay/replay_lme_planning.csv")
Q <- read.csv("D:/2020_RiskyReplay/results/questionnaire_results.csv")

#######################################
# CLEAN DATA
#######################################

cleandata <- function(d) {

  d$Subject <- as.factor(d$Subject)
  d$Acc <- as.factor(d$Acc)
  
  d$P <- round(d$P,1)
  d$Certainty <- rep(NA,nrow(d))
  d$Certainty[d$P==.5] <- 1
  d$Certainty[d$P==.3 | d$P==.7] <- 2
  d$Certainty[d$P==.1 | d$P==.9] <- 3
  
  d$RewProb <- rep(NA,nrow(d))
  d$RewProb[d$nV_1>d$nV_2] <- round(d$P[d$nV_1>d$nV_2],1)
  d$RewProb[d$nV_1<d$nV_2] <- round(1-d$P[d$nV_1<d$nV_2],1)
  d$RewProb[d$nV_1>0 & d$nV_2>0] <- 1
  d$RewProb[d$nV_1<0 & d$nV_2<0] <- 0
  
  d$RewVal <- rep(NA,nrow(d))
  d$RewVal[d$nV_1>d$nV_2] <- d$nV_1[d$nV_1>d$nV_2]
  d$RewVal[d$nV_1<d$nV_2] <- d$nV_2[d$nV_1<d$nV_2]
  
  d$LossVal <- rep(NA,nrow(d))
  d$LossVal[d$nV_1<d$nV_2] <- d$nV_1[d$nV_1<d$nV_2]
  d$LossVal[d$nV_1>d$nV_2] <- d$nV_2[d$nV_1>d$nV_2]
  
  if (any(names(d)=="Replay_rewarding")){
    d$Replay_average <- (d$Replay_rewarding + d$Replay_aversive)/2
  }
  
  if (!any(names(d)=="OverallBlock")){
    d$OverallBlock <- rep(NA,nrow(d))
    d$OverallBlock[d$Practice==1] <- 0
    d$OverallBlock[d$Practice==0] <- d$Block[d$Practice==0]
  }
  
  d$ExpHalf <- rep(NA,nrow(d))
  d$ExpHalf[d$OverallBlock <= 5] <- "first"
  d$ExpHalf[d$OverallBlock > 5] <- "second"
  
  d <- d %>%
    group_by(Subject,OverallBlock) %>%
    mutate(BlockAcc=mean(Acc=="correct")) %>%
    group_by(Subject) %>%
    mutate(SubjAcc=mean(Acc=="correct"))
  
  d <- d %>% filter(RewProb!=0,RewProb!=1)
  
  # add in cumulative path experience
  d$PathTransition <- rep("safe",nrow(d))
  d$PathTransition[d$Choice==1 & d$Outcome>0] <- "rewarding"
  d$PathTransition[d$Choice==1 & d$Outcome<0] <- "aversive"
  
  d <- d %>%
    group_by(Subject,Lag) %>%
    mutate(RewPathExp_overall = cumsum(PathTransition=="rewarding"),
           LossPathExp_overall = cumsum(PathTransition=="aversive")) %>%
    ungroup() %>%
    group_by(Subject,Lag,ExpHalf) %>%
    mutate(RewPathExp_half = cumsum(PathTransition=="rewarding"),
           LossPathExp_half = cumsum(PathTransition=="aversive")) %>%
    ungroup() %>%
    group_by(Subject,Lag,OverallBlock) %>%
    mutate(RewPathExp_block = cumsum(PathTransition=="rewarding"),
           LossPathExp_block = cumsum(PathTransition=="aversive"))
  
  # Shift replay to NEXT trial
  d <- d %>%
    group_by(Subject,Lag,OverallBlock) %>%
    mutate(PrevReplay_differential = c(NA,Replay_differential[1:(n()-1)]),
           PrevReplay_rewarding = c(NA,Replay_rewarding[1:(n()-1)]),
           PrevReplay_aversive = c(NA,Replay_aversive[1:(n()-1)]))
  
  # finish
  return(as.data.frame(d))
}

d <- cleandata(d)

#######################################
# PCA ON QUESTIONNAIRES
#######################################

Q <- Q[Q$Subject %in% intersect(unique(d$Subject),Q$Subject),]

pca <- prcomp(Q[,c(2,3,5,6,7,8,9)],center=TRUE,scale.=TRUE)
eig <- summary(pca)
PC <- as.data.frame(pca$x)

# # scree plot
# screeplot(pca,type="l")
# abline(h=1,col="red",lty=5)
# 
# # correlation matrix
# library(ggcorrplot)
# C <- cbind(Q,PC) %>% select(starts_with("score_"),starts_with("PC",ignore.case=FALSE))
# corr <- round(cor(C),1)
# p.mat <- cor_pmat(C)
# ggcorrplot(corr,type="upper",lab=FALSE)


addQ <- function(d,Q){
  
  d$IUS <- rep(NA,nrow(d))
  d$Worry <- rep(NA,nrow(d))
  d$Risk <- rep(NA,nrow(d))
  d$PCA1 <- rep(NA,nrow(d))
  d$PCA2 <- rep(NA,nrow(d))
  for (subject in Q$Subject) {
    idx <- d$Subject==subject
    d$IUS[idx] <- Q$pcnt_IUS[Q$Subject==subject]
    d$Worry[idx] <- Q$pcnt_worry[Q$Subject==subject]
    d$Risk[idx] <- Q$pcnt_risk[Q$Subject==subject]
    d$PCA1[idx] <- PC$PC1[Q$Subject==subject]
    d$PCA2[idx] <- PC$PC2[Q$Subject==subject]
  }
  
  d$Anxiety <- rowMeans(select(d,IUS,Worry))
  return(d)
}

d <- addQ(d,Q)

# # Plot correlations
# pd <- cbind(Q,PC) %>%
#   gather(PCAnum,PCAval,PC1:PC2)
#   
# p1 <- ggplot(data=pd,aes(x=score_IUS,y=PCAval,group=PCAnum)) + 
#   geom_point() + 
#   geom_smooth(method=lm) + 
#   facet_wrap(~PCAnum) + theme_classic() + 
#   coord_cartesian(xlim=c(10,55)) + scale_x_continuous(breaks=seq(10,55,5))
# 
# ggplot(data=pd,aes(x=score_worry,y=PCAval,group=PCAnum)) + 
#   geom_point() + 
#   geom_smooth(method=lm) + 
#   facet_wrap(~PCAnum) + theme_classic()
# 
# ggplot(data=pd,aes(x=score_risk,y=PCAval,group=PCAnum)) + 
#   geom_point() + 
#   geom_smooth(method=lm) + 
#   facet_wrap(~PCAnum) + theme_classic()

#######################################
# BEHAVIOURAL ANALYSIS
#######################################

md <- filter(d,Lag==d$Lag[1])
md$EV <- scale(md$EV,center=TRUE,scale=FALSE)
md$PCA1 <- scale(md$PCA1,center=TRUE,scale=FALSE)
md$PCA2 <- scale(md$PCA2,center=TRUE,scale=FALSE)
md$RT <- scale(md$RT,center=TRUE,scale=FALSE)
md$Trial <- scale(md$Trial,center=TRUE,scale=FALSE)
md$RewVal <- scale(md$RewVal,center=TRUE,scale=FALSE)
md$LossVal <- scale(md$LossVal,center=TRUE,scale=FALSE)
md$RewProb <- scale(md$RewProb,center=TRUE,scale=FALSE)

md$RewPathExp_block <- scale(md$RewPathExp_block,center=TRUE,scale=FALSE)
md$LossPathExp_block <- scale(md$LossPathExp_block,center=TRUE,scale=FALSE)
md$DiffPathExp <- scale(md$RewPathExp_block - md$LossPathExp_block,center=TRUE,scale=FALSE)

m0 <- glmer(Choice ~ EV + Certainty + RT + (1|Subject),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m1 <- glmer(Choice ~ RewVal*LossVal*RewProb + Certainty + RT + (1|Subject),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

summary(m1)




m2 <- glmer(Choice ~ RewPathExp_block*LossPathExp_block + RewVal*LossVal*RewProb + Certainty + RT + (1|Subject),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))
summary(m2)

m3 <- glmer(Choice ~ DiffPathExp + RewVal*LossVal*RewProb + Certainty + RT + (1|Subject),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))
summary(m3)


p1 <- interact_plot(m1,pred=RewProb,modx=RewVal,modx.values=quantile(md$RewVal,c(.05,.5,.95)),
              interval=TRUE,vary.lty=FALSE) + 
  theme_classic()

p2 <- interact_plot(m1,pred=RewProb,modx=LossVal,modx.values=quantile(md$LossVal,c(.05,.5,.95)),
                    interval=TRUE,vary.lty=FALSE) + 
  theme_classic()

g <- grid.arrange(p1,p2,nrow=1)

ggsave("rewvalrewprob.svg",g)


summary(m0)
pred <- jtools::make_predictions(m0,pred="EV")

ggplot(data=pred, aes(x=EV,y=Choice)) + 
  geom_ribbon(aes(ymax=ymax,ymin=ymin),alpha=.2) + 
  geom_line(size=1) +
  theme_classic() +
  coord_cartesian(xlim=c(-9,11)) + scale_x_continuous(breaks=seq(-9,11,2))

m1 <- glmer(Choice ~ EV + PCA1 + Certainty + (1|Subject),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m2 <- glmer(Choice ~ EV + PCA2 + Certainty + (1|Subject),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m3 <- glmer(Choice ~ EV + PCA1 + PCA2 + Certainty + (1|Subject),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m4 <- glmer(Choice ~ EV*PCA1 + Certainty + (1|Subject),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m5 <- glmer(Choice ~ EV*PCA2 + Certainty + (1|Subject),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m6 <- glmer(Choice ~ EV*PCA1 + EV*PCA2 + Certainty + (1|Subject),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m7 <- glmer(Choice ~ EV + Certainty*PCA1 + (1|Subject),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m8 <- glmer(Choice ~ EV + Certainty*PCA2 + (1|Subject),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m9 <- glmer(Choice ~ EV + Certainty*PCA1 + Certainty*PCA2 + (1|Subject),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m10 <- glmer(Choice ~ EV*PCA1 + Certainty*PCA1 + (1|Subject),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m11 <- glmer(Choice ~ EV*PCA2 + Certainty*PCA2 + (1|Subject),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m12 <- glmer(Choice ~ EV*PCA1 + Certainty*PCA1 + EV*PCA2 + Certainty*PCA2 + (1|Subject),
             data=md, family="binomial",
             control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

anova(m0,m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12)

#######################################
# REPLAY ANALYSIS (VALUE)
#######################################

# Rearrange so that replay is spread out across rows instead of columns
md <- d %>%
  gather(Replay_type, Sequenceness, c(Replay_rewarding:Replay_differential,Replay_average), factor_key=TRUE) %>%
  mutate(Replay_type=gsub("Replay_","",Replay_type))
md$Replay_type <- as.factor(md$Replay_type)

md$PathVal <- rep(NA,nrow(md))
md$PathVal[md$Replay_type=="rewarding"] <- apply(select(filter(md,Replay_type=="rewarding"),nV_1:nV_2),1,max)
md$PathVal[md$Replay_type=="aversive"] <- apply(select(filter(md,Replay_type=="aversive"),nV_1:nV_2),1,min)

md$PathProb <- rep(NA,nrow(md))
md$PathProb[md$Replay_type=="rewarding"] <- md$RewProb[md$Replay_type=="rewarding"]
md$PathProb[md$Replay_type=="aversive"] <- round(1-md$RewProb[md$Replay_type=="aversive"],1)

md$PathCertainty <- md$PathProb
md$PathCertainty[which(md$PathProb==0.1)] <- rep(0.9,sum(md$PathProb==0.1,na.rm=TRUE))
md$PathCertainty[which(md$PathProb==0.3)] <- rep(0.7,sum(md$PathProb==0.1,na.rm=TRUE))

md$PathNum <- rep(NA,nrow(md))
md$PathNum[md$Replay_type=="rewarding" & md$nV_1>md$nV_2] <- 1
md$PathNum[md$Replay_type=="rewarding" & md$nV_1<md$nV_2] <- 2
md$PathNum[md$Replay_type=="aversive" & md$nV_1<md$nV_2] <- 1
md$PathNum[md$Replay_type=="aversive" & md$nV_1>md$nV_2] <- 2
md$PathNum <- as.factor(md$PathNum)

md$PathExp <- rep(NA,nrow(md))
md$PathExp[md$Replay_type=="rewarding"] <- md$RewPathExp_block[md$Replay_type=="rewarding"]
md$PathExp[md$Replay_type=="aversive"] <- md$LossPathExp_block[md$Replay_type=="aversive"]
md$PathExp <- log(md$PathExp+1)

md$PathRecency <- rep(NA,nrow(md))
md$PathRecency[md$Replay_type=="rewarding"] <- md$RewPathRecency_block[md$Replay_type=="rewarding"]
md$PathRecency[md$Replay_type=="aversive"] <- md$LossPathRecency_block[md$Replay_type=="aversive"]
md$PathRecency <- log(md$PathRecency+1)

md$PrevPathType <- md$PrevOutcome > 0

md$SelectedPath <- (md$PathTransition=="rewarding" & md$Replay_type=="rewarding") | (md$PathTransition=="aversive" & md$Replay_type=="aversive")

md$PrevSelected <- rep(NA,nrow(md))
md$PrevSelected[md$PrevOutcome>0 & md$Replay_type=="rewarding"] = TRUE
md$PrevSelected[md$PrevOutcome>0 & md$Replay_type=="aversive"] = FALSE
md$PrevSelected[md$PrevOutcome<0 & md$Replay_type=="rewarding"] = FALSE
md$PrevSelected[md$PrevOutcome<0 & md$Replay_type=="aversive"] = TRUE
md$PrevSelected <- as.factor(md$PrevSelected)

md$DiffPathExp <- md$RewPathExp_block - md$LossPathExp_block

md$Choice <- as.factor(md$Choice)

md$UPE <- abs(md$PEbyEV) # unsigned prediction error

md <- filter(md,
             Lag>10,Lag<100,
             Replay_type=="rewarding" | Replay_type=="aversive")

# md <- md %>%
#   group_by(Subject,ExpTrial) %>%
#   mutate(AvSequenceness=mean(Sequenceness)) %>%
#   ungroup() %>%
#   mutate(Outlier=abs(scale(AvSequenceness))>3)
# md <- filter(md,!Outlier)

md$PathExp <- scale(md$PathExp,center=TRUE,scale=FALSE)
md$DiffPathExp <- scale(md$DiffPathExp,center=TRUE,scale=FALSE)
md$PathVal <- scale(md$PathVal,center=TRUE,scale=FALSE)
md$RT <- scale(md$RT,center=TRUE,scale=FALSE)
md$PathProb <- scale(md$PathProb,center=TRUE,scale=FALSE)
md$PathVal <- scale(md$PathVal,center=TRUE,scale=FALSE)
md$PathRecency <- scale(md$PathRecency,center=TRUE,scale=FALSE)

md$CatPEbyEV <- md$PEbyEV > 0

md$PEbyEV <- scale(md$PEbyEV,center=TRUE,scale=FALSE)
md$Outcome <- scale(md$Outcome,center=TRUE,scale=FALSE)
md$UPE <- scale(md$UPE,center=TRUE,scale=FALSE)

md$PrevOutcome <- scale(md$PrevOutcome,center=TRUE,scale=FALSE)
md$PEbyRew <- scale(md$PEbyRew,center=TRUE,scale=FALSE)
md$PEbyLoss <- scale(md$PEbyLoss,center=TRUE,scale=FALSE)



m0 <- lmer(Sequenceness ~ Replay_type + RT + (1|Subject/Lag),
            data=md, REML=FALSE,
            control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))
summary(m0)

m1 <- lmer(Sequenceness ~ Replay_type*Choice + RT + (1|Subject/Lag),
           data=md, REML=FALSE,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))
summary(m1)

m2 <- lmer(Sequenceness ~ PathExp*Replay_type*Choice + RT + (1|Subject/Lag),
           data=md, REML=FALSE,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))
summary(m2)
# 
# pred <- jtools::make_predictions(m2,pred="PathExp",
#                            at=list(Replay_type=unique(md$Replay_type),
#                                                  Choice=unique(md$Choice))) %>%
#         mutate(semupper=Sequenceness+((ymax-ymin)/3.96),
#                semlower=Sequenceness-((ymax-ymin)/3.96))
# 
# ggplot(pred,aes(x=PathExp,y=Sequenceness,group=Replay_type,color=Replay_type,fill=Replay_type)) + 
#   geom_ribbon(aes(ymin=semlower,ymax=semupper),alpha=.2,color=NA) + 
#   geom_line(size=1) + 
#   facet_wrap(~Choice) + 
#   theme_classic() + 
#   coord_cartesian(xlim=c(-1,1.5),ylim=c(-0,0.045)) + 
#   scale_y_continuous(breaks=seq(-0.005,0.045,0.005)) + 
#   scale_fill_manual(values=c("#FF0074","#00FF9B")) + 
#   scale_color_manual(values=c("#FF0074","#00FF9B"))

# effect of path probability
mA <- lmer(Sequenceness ~ PathProb*Choice + RT + (1|Subject/Lag),
           data=md, REML=FALSE,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))
summary(mA)

pred <- jtools::make_predictions(mA,pred="PathProb",
                                 at=list(Choice=unique(md$Choice))) %>%
  mutate(semupper=Sequenceness+((ymax-ymin)/3.96),
         semlower=Sequenceness-((ymax-ymin)/3.96))

ggplot(pred,aes(x=PathProb,y=Sequenceness,group=Choice,color=Choice,fill=Choice)) + 
  geom_ribbon(aes(ymin=semlower,ymax=semupper),alpha=.2,color=NA) + 
  geom_line(size=1) + 
  theme_classic() + 
  coord_cartesian(ylim=c(-0,0.04))

# effect of path recency (higher numbers = experienced longer ago)
m0 <- lmer(Sequenceness ~ PathRecency*Replay_type*PathProb + RT + (1|Subject/Lag),
           data=md, REML=TRUE,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))
summary(m0)

pred <- jtools::make_predictions(m0,pred="PathRecency",
                                 at=list(Replay_type=unique(md$Replay_type),
                                         PathProb=unique(md$PathProb))) %>%
  mutate(semupper=Sequenceness+((ymax-ymin)/3.96),
         semlower=Sequenceness-((ymax-ymin)/3.96))

ggplot(pred,aes(x=PathRecency,y=Sequenceness,
                group=Replay_type,color=PathProb)) + 
  geom_point() + 
  theme_classic() + 
  coord_cartesian(ylim=c(0,0.04))

# compare reward & loss paths
m0 <- lmer(Sequenceness ~ Replay_type + RT + (1|Subject/Lag),
           data=md, REML=TRUE,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))
summary(m0)

m1 <- lmer(Sequenceness ~ Replay_type*Choice + RT + (1|Subject/Lag),
           data=md, REML=TRUE,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))
summary(m1)

cat_plot(m1,pred=Choice,modx=Replay_type)
emmeans(m1,pairwise~Choice*Replay_type)


# compare to rescorla-wagner learned value
alpha = 5  # 1, 3, 5, 7, or 9

md$RS_pathval <- rep(NA,nrow(md))
md$RS_pathval[md$PathNum==1] <- select(md[md$PathNum==1,],contains(paste("RW_a",alpha,"_path1",sep="")))[,1]
md$RS_pathval[md$PathNum==2] <- select(md[md$PathNum==2,],contains(paste("RW_a",alpha,"_path2",sep="")))[,1]

md$RS_replay_type <- rep(NA,nrow(md))
md$RS_replay_type[md$RS_pathval>0] <- "rewarding"
md$RS_replay_type[md$RS_pathval<1] <- "aversive"

m1 <- lmer(Sequenceness ~ PathVal*Choice + RT + (1|Subject/Lag),
           data=md, REML=TRUE,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))
summary(m1)

m2 <- lmer(Sequenceness ~ RS_replay_type*Choice + RT + (1|Subject/Lag),
           data=md, REML=TRUE,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))
summary(m2)

anova(m1,m2)


# effect of previous trial prediction error
m0 <- lmer(Sequenceness ~ PEbyEV*Replay_type + RT + (1|Subject/Lag),
           data=md, REML=TRUE,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))
summary(m0)

interact_plot(m0,pred=PEbyEV,modx=Replay_type,interval=TRUE,colors=c("#FF0074","#00FF9B"),vary.lty=FALSE) + theme_classic()

m1 <- lmer(Sequenceness ~ PEbyEV_e1*Replay_type*Choice + RT + (1|Subject/Lag),
           data=md, REML=TRUE,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m2 <- lmer(Sequenceness ~ PEbyEV_e2*Replay_type*Choice + RT + (1|Subject/Lag),
           data=md, REML=TRUE,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m3 <- lmer(Sequenceness ~ PEbyEV_e3*Replay_type*Choice + RT + (1|Subject/Lag),
           data=md, REML=TRUE,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m4 <- lmer(Sequenceness ~ PEbyEV_e4*Replay_type*Choice + RT + (1|Subject/Lag),
           data=md, REML=TRUE,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m5 <- lmer(Sequenceness ~ PEbyEV_e5*Replay_type*Choice + RT + (1|Subject/Lag),
           data=md, REML=TRUE,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m6 <- lmer(Sequenceness ~ PEbyEV_e6*Replay_type*Choice + RT + (1|Subject/Lag),
           data=md, REML=TRUE,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m7 <- lmer(Sequenceness ~ PEbyEV_e7*Replay_type*Choice + RT + (1|Subject/Lag),
           data=md, REML=TRUE,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m8 <- lmer(Sequenceness ~ PEbyEV_e8*Replay_type*Choice + RT + (1|Subject/Lag),
           data=md, REML=TRUE,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m9 <- lmer(Sequenceness ~ PEbyEV_e9*Replay_type*Choice + RT + (1|Subject/Lag),
           data=md, REML=TRUE,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m10 <- lmer(Sequenceness ~ PEbyEV_e10*Replay_type*Choice + RT + (1|Subject/Lag),
           data=md, REML=TRUE,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m11 <- lmer(Sequenceness ~ PEbyEV_e11*Replay_type*Choice + RT + (1|Subject/Lag),
           data=md, REML=TRUE,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m12 <- lmer(Sequenceness ~ PEbyEV_e12*Replay_type*Choice + RT + (1|Subject/Lag),
           data=md, REML=TRUE,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))
#######################################
# REPLAY ANALYSIS (CHOICE)
#######################################

# First check that all subjects have an appropriate EV curve
for (subject in unique(d$Subject)) { 
  
  md <- d %>%
    filter(Lag>10, Lag<100, Subject==subject)
  
  md$EV <- scale(md$EV,center=TRUE,scale=FALSE)
  
  m <- glmer(Choice ~ EV + (1|Lag),
             data=md, family="binomial",
             control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))
  
  pred <- jtools::make_predictions(m,pred="EV")
  pred$Subject <- rep(subject,nrow(pred))
  pred$PCA1 <- rep(unique(md$PCA1[md$Subject==subject]),nrow(pred))
  pred$PCA2 <- rep(unique(md$PCA2[md$Subject==subject]),nrow(pred))
  
  pred <- select(pred,Choice,EV,ymax,ymin,Subject,PCA1,PCA2)
  
  m.summary <- summary(m)
  thisbeta <- data.frame(Subject=subject,
                         Beta=fixef(m)[2],
                         Zval=m.summary$coefficients[2,3],
                         Pval=m.summary$coefficients[2,4])
  
  if (subject==unique(d$Subject)[1]) {
    allpred <- pred
    allbetas <- thisbeta
  } else {
    allpred <- rbind(allpred,pred)
    allbetas <- rbind(allbetas,thisbeta)
  }
  
}

ggplot(data=allpred, aes(x=EV,y=Choice,group=Subject,color=Subject)) + 
  geom_line(size=1) + 
  theme_classic() + 
  scale_colour_viridis(discrete=TRUE)





md <- d %>%
  filter(Lag>10, Lag<100) # Subject!=97403

# md <- md %>%
#   group_by(Subject) %>%
#   mutate(zReplay = abs(scale(Replay_differential))) %>%
#   group_by(Subject,ExpTrial) %>%
#   mutate(outliers = any(zReplay>3)) %>%
#   filter(outliers==FALSE)

md$BinDiffReplay <- md$Replay_differential > 0
md$BinDiffReplay[md$Replay_differential>0] <- rep("rewarding",sum(md$Replay_differential>0))
md$BinDiffReplay[md$Replay_differential<0] <- rep("aversive",sum(md$Replay_differential<0))

# md <- md %>%
#   group_by(Subject,ExpTrial) %>%
#   mutate(AvSequenceness=mean(Replay_differential)) %>%
#   ungroup() %>%
#   mutate(Outlier=abs(scale(AvSequenceness))>3)
# md <- filter(md,!Outlier)

md$origEV <- md$EV
md$EV <- scale(md$EV,center=TRUE,scale=FALSE)
md$RT <- scale(md$RT,center=TRUE,scale=FALSE)
# md$Replay_differential <- scale(md$Replay_differential,center=TRUE,scale=FALSE)
# md$Replay_average <- scale(md$Replay_average,center=TRUE,scale=FALSE)
md$Replay_rewarding <- scale(md$Replay_rewarding,center=TRUE,scale=FALSE)
md$Replay_aversive <- scale(md$Replay_aversive,center=TRUE,scale=FALSE)
md$RewVal <- scale(md$RewVal,center=TRUE,scale=FALSE)
md$LossVal <- scale(abs(md$LossVal),center=TRUE,scale=FALSE)
md$RewProb <- scale(md$RewProb,center=TRUE,scale=FALSE)
md$DiffPathExp <- scale(md$RewPathExp_block - md$LossPathExp_block,center=TRUE,scale=FALSE)
md$PEbyEV <- scale(md$PEbyEV,center=TRUE,scale=FALSE)
md$PrevOutcome <- scale(md$PrevOutcome,center=TRUE,scale=FALSE)

# lm0 <- lmer(Replay_differential ~ Choice + RT + Certainty + Replay_average + (1|Subject/Lag),
#             data=md, REML=TRUE,
#             control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))
# 
# lm1 <- lmer(Replay_differential ~ RewVal*LossVal*RewProb + RT + Certainty + Replay_average + (1|Subject/Lag),
#             data=md, REML=TRUE,
#             control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))
# 
# summary(lm1)
# 
# 
# g <- interact_plot(lm1,pred=RewVal,modx=LossVal,modx.values=quantile(md$LossVal,c(.25,.5,.75)),
#               interval=TRUE,vary.lty=FALSE,colors="Reds") + 
#   theme_classic() + 
#   coord_cartesian(ylim=c(-0.02,0.03)) + 
#   scale_y_continuous(breaks=seq(-0.02,0.03,0.01))
# 
# ggsave("rewval_replay.svg",g)
# 
# 
# g <- jtools::effect_plot(lm1,pred=RewProb,interval=TRUE) + 
#   theme_classic() + 
#   coord_cartesian(ylim=c(-0.02,0.03)) + 
#   scale_y_continuous(breaks=seq(-0.02,0.03,0.01))
# ggsave("rewprob_replay.svg",g)





m0 <- glmer(Choice ~ EV + Certainty + RT + (1|Subject/Lag),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m1 <- glmer(Choice ~ Replay_differential + RT + (1|Subject/Lag),
           data=md, family="binomial",
           control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m2 <- glmer(Choice ~ EV*Replay_differential + Certainty + RT + (1|Subject/Lag),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))




# (prediction error analysis)
# d$PEbyEV[md$PEbyEV=="0"] = NA # ignore trials where the previous outcome was the avoidance +1 point

m2a <- glmer(Choice ~ EV*Replay_differential*PEbyEV + Certainty + RT + (1|Subject/Lag),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m2b <- glmer(Choice ~ EV*Replay_differential*PrevOutcome + Certainty + RT + (1|Subject/Lag),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m2c <- glmer(Choice ~ EV*Replay_differential*PEbyEV + EV*Replay_differential*PrevOutcome + Certainty + RT + (1|Subject/Lag),
             data=md, family="binomial",
             control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

interact_plot(m2a,pred=PEbyEV,modx=Replay_differential,modx.values=quantile(md$Replay_differential,c(.05,.95)),interval=TRUE,
              vary.lty=FALSE,colors=c("#FF006C","#00FF7C")) + theme_classic() + coord_cartesian(ylim=c(0,1))

interact_plot(m2b,pred=PrevOutcome,modx=Replay_differential,modx.values=quantile(md$Replay_differential,c(.05,.95)),interval=TRUE,
              vary.lty=FALSE,colors=c("#FF006C","#00FF7C")) + theme_classic() + coord_cartesian(ylim=c(0,1))




m3 <- glmer(Choice ~ EV*Replay_rewarding + Certainty + RT + (1|Subject/Lag),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m4 <- glmer(Choice ~ EV*Replay_aversive + Certainty + RT + (1|Subject/Lag),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m5 <- glmer(Choice ~ EV*Replay_rewarding+ EV*Replay_aversive + Certainty + RT + (1|Subject/Lag),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

# m6 <- glmer(Choice ~ EV*Replay_rewarding*Replay_aversive + Certainty + RT + (1|Subject/Lag),
#             data=md, family="binomial",
#             control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

anova(m0,m1)
summary(m1)


# (see how replay & EV differ, as function of choice)
md <- d %>%
  filter(Lag>10, Lag<100, Choice==0)

md$EVacc <- md$EV
md$EVacc[md$EV>1 & md$Choice==0] <- md$EV[md$EV>1 & md$Choice==0] * (-1)
md$EVacc[md$EV<1 & md$Choice==0] <- md$EV[md$EV<1 & md$Choice==0] * (-1)

md$EVacc <- scale(md$EVacc,center=TRUE,scale=FALSE)
md$RT <- scale(md$RT,center=TRUE,scale=FALSE)
md$Certainty <- scale(md$Certainty,center=TRUE,scale=FALSE)

m6 <- lmer(Replay_rewarding ~ EV + Certainty + RT + (1|Subject/Lag),
            data=md, REML=FALSE,
            control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))
pd <- md %>% 
  mutate(rEV = round(EV)) %>% 
  group_by(Subject,rEV) %>% 
  mutate(N=n(),Replay=mean(Replay_rewarding)) %>% 
  group_by(rEV) %>% 
  summarise(Replay=mean(Replay),N=mean(N))
jtools::effect_plot(m6,pred=EV,interval=TRUE) + 
  geom_point(data=pd,aes(x=rEV,y=Replay,size=N))

m7 <- lmer(Replay_aversive ~ EV + Certainty + RT + (1|Subject/Lag),
           data=md, REML=FALSE,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))
pd <- md %>% 
  mutate(rEV = round(EV)) %>% 
  group_by(Subject,rEV) %>% 
  mutate(N=n(),Replay=mean(Replay_aversive)) %>% 
  group_by(rEV) %>% 
  summarise(Replay=mean(Replay),N=mean(N))
jtools::effect_plot(m7,pred=EV,interval=TRUE) + 
  geom_point(data=pd,aes(x=rEV,y=Replay,size=N))



md <- d %>%
  filter(Lag>10, Lag<100, Choice==1)
md$EV <- scale(md$EV,center=TRUE,scale=FALSE)
md$RT <- scale(md$RT,center=TRUE,scale=FALSE)
md$Certainty <- scale(md$Certainty,center=TRUE,scale=FALSE)

m8 <- lmer(Replay_rewarding ~ EV + Certainty + RT + (1|Subject/Lag),
           data=md, REML=FALSE,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))
pd <- md %>% 
  mutate(rEV = round(EV)) %>% 
  group_by(Subject,rEV) %>% 
  mutate(N=n(),Replay=mean(Replay_rewarding)) %>% 
  group_by(rEV) %>% 
  summarise(Replay=mean(Replay),N=mean(N))
jtools::effect_plot(m8,pred=EV,interval=TRUE) + 
  geom_point(data=pd,aes(x=rEV,y=Replay,size=N))

m9 <- lmer(Replay_aversive ~ EV + Certainty + RT + (1|Subject/Lag),
           data=md, REML=FALSE,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))
pd <- md %>% 
  mutate(rEV = round(EV)) %>% 
  group_by(Subject,rEV) %>% 
  mutate(N=n(),Replay=mean(Replay_aversive)) %>% 
  group_by(rEV) %>% 
  summarise(Replay=mean(Replay),N=mean(N))
jtools::effect_plot(m9,pred=EV,interval=TRUE) + 
  geom_point(data=pd,aes(x=rEV,y=Replay,size=N))






md <- d %>%
  filter(Lag>10, Lag<100, Choice==1)
md$EV <- scale(md$EV,center=TRUE,scale=FALSE)
md$RT <- scale(md$RT,center=TRUE,scale=FALSE)
md$Certainty <- scale(md$Certainty,center=TRUE,scale=FALSE)

m7 <- lmer(Replay_rewarding ~ EV + Certainty + RT + (1|Subject/Lag),
           data=md, REML=FALSE,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))
pd <- md %>% 
  mutate(rEV = round(EV)) %>% 
  group_by(Subject,rEV) %>% 
  mutate(N=n(),Replay=mean(Replay_rewarding)) %>% 
  group_by(rEV) %>% 
  summarise(Replay=mean(Replay),N=mean(N))
jtools::effect_plot(m7,pred=EV,interval=TRUE) + 
  geom_point(data=pd,aes(x=rEV,y=Replay,size=N))




# get indifference points
pred <- jtools::make_predictions(m2,pred="EV")
pred$EV[abs(pred$Choice - 0.5) == min(abs(pred$Choice - 0.5))]

pred <- jtools::make_predictions(m2,pred="EV",
        at=list(Replay_differential=quantile(md$Replay_differential,c(.05,.95))))
idx1 <- pred$Replay_differential==min(pred$Replay_differential)
idx2 <- pred$Replay_differential==max(pred$Replay_differential)
pred$Replay_differential[idx1] <- rep("aversive",sum(idx1))
pred$Replay_differential[idx2] <- rep("rewarding",sum(idx2))

tmp <- pred[idx1,]
tmp$EV[abs(tmp$Choice - 0.5) == min(abs(tmp$Choice - 0.5))]

tmp <- pred[idx2,]
tmp$EV[abs(tmp$Choice - 0.5) == min(abs(tmp$Choice - 0.5))]
                                                                     
# estimated means
pred <- jtools::make_predictions(m2,pred="Replay_differential")
mean(pred$Replay_differential[pred$Choice>0.5])



                                                                     
acf(resid(m1))
durbinWatsonTest(resid(m1)) # values close to 2 indicate no autocorrelation (<2 = positive, >2 = negative) - rule of thumb is 1.5 to 2.5 is fine



# evplots <- function(pred) {
#   
#   pred_choice <- pred %>%
#     select(Choice,Replay_differential,EV) %>%
#     spread(Replay_differential,Choice)
#   names(pred_choice) <- c("EV","Aversive_Choice","Rewarding_Choice")
#   
#   pred_ymin <- pred %>%
#     select(ymin,Replay_differential,EV) %>%
#     spread(Replay_differential,ymin)
#   names(pred_ymin) <- c("EV","Aversive_ymin","Rewarding_ymin")
#   
#   pred_ymax <- pred %>%
#     select(ymax,Replay_differential,EV) %>%
#     spread(Replay_differential,ymax)
#   names(pred_ymax) <- c("EV","Aversive_ymax","Rewarding_ymax")
#   
#   pred <- cbind(pred_choice,pred_ymin[,2:3],pred_ymax[,2:3])
#   
#   pred$Aversive_yminse <- pred$Aversive_Choice - (pred$Aversive_ymax-pred$Aversive_ymin)/3.96
#   pred$Aversive_ymaxse <- pred$Aversive_Choice + (pred$Aversive_ymax-pred$Aversive_ymin)/3.96
#   
#   pred$Rewarding_yminse <- pred$Rewarding_Choice - (pred$Rewarding_ymax-pred$Rewarding_ymin)/3.96
#   pred$Rewarding_ymaxse <- pred$Rewarding_Choice + (pred$Rewarding_ymax-pred$Rewarding_ymin)/3.96
#   
#   
#   return(pred)
# }

pred <- jtools::make_predictions(m2,
                                 pred="EV",modx="Replay_differential",
                                 at=list(Replay_differential=quantile(d$Replay_differential,c(.05,.5,.95))))
names(pred)[names(pred)=="Replay_differential"] <- "Replay"
write.csv(pred,"ev_replay_differential.csv")

pred <- jtools::make_predictions(m5,
                                 pred="EV",modx="Replay_rewarding",
                                 at=list(Replay_rewarding=quantile(d$Replay_rewarding,c(.05,.5,.95))))
names(pred)[names(pred)=="Replay_rewarding"] <- "Replay"
write.csv(pred,"ev_replay_reward.csv")

pred <- jtools::make_predictions(m5,
                                 pred="EV",modx="Replay_aversive",
                                 at=list(Replay_aversive=quantile(d$Replay_aversive,c(.05,.5,.95))))
names(pred)[names(pred)=="Replay_aversive"] <- "Replay"
write.csv(pred,"ev_replay_loss.csv")


interact_plot(m1,data=md,pred=EV,modx=Replay_differential,modx.values=quantile(d$Replay_differential,c(.05,.5,.95)),
              interval=TRUE,colors=c("#FF0055","#737373","#00FF8F"),vary.lty=FALSE) + 
  theme_classic()

ss <- sim_slopes(m1,pred=EV,modx=Replay_differential,modx.values=quantile(d$Replay_differential,c(.05,.95)),
                 cond.int=TRUE)


pd <- md %>%
  mutate(Choice=as.factor(Choice)) %>%
  group_by(Subject,Choice) %>%
  summarise(Replay=mean(Replay_differential))

ggplot(data=pd,aes(x=Choice,y=Replay,group=Subject,color=Choice)) + 
  geom_line(alpha=.5,colour="black") + 
  geom_point(size=3) + 
  theme_classic()

write.csv(pd,"choice_replay.csv")


m1 <- glmer(Choice ~ EV*Replay_differential + Certainty + RT + (1|Subject/Lag),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m2 <- glmer(Choice ~ EV*Replay_differential + PCA1 + Certainty + (1|Subject/Lag),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m3 <- glmer(Choice ~ EV*Replay_differential + PCA2 + Certainty + (1|Subject/Lag),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m4 <- glmer(Choice ~ EV*Replay_differential + PCA1 + PCA2 + Certainty + (1|Subject/Lag),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m5 <- glmer(Choice ~ EV*Replay_differential*PCA1 + Certainty + (1|Subject/Lag),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m6 <- glmer(Choice ~ EV*Replay_differential*PCA1 + PCA2 + Certainty + (1|Subject/Lag),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m7 <- glmer(Choice ~ EV*Replay_differential*PCA2 + Certainty + (1|Subject/Lag),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m8 <- glmer(Choice ~ EV*Replay_differential*PCA2 + PCA1 + Certainty + (1|Subject/Lag),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m9 <- glmer(Choice ~ EV*Replay_differential*PCA1 + EV*Replay_differential*PCA2 + Certainty + (1|Subject/Lag),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m10 <- glmer(Choice ~ EV*Replay_rewarding*PCA1 + EV*Replay_rewarding*PCA2 + Certainty + (1|Subject/Lag),
             data=md, family="binomial",
             control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m11 <- glmer(Choice ~ EV*Replay_aversive*PCA1 + EV*Replay_aversive*PCA2 + Certainty + (1|Subject/Lag),
             data=md, family="binomial",
             control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m12 <- glmer(Choice ~ EV*Replay_rewarding*PCA1 + EV*Replay_rewarding*PCA2 + EV*Replay_aversive*PCA1 + EV*Replay_aversive*PCA2 + Certainty + (1|Subject/Lag),
             data=md, family="binomial",
             control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

anova(m1,m2,m3,m4,m5,m6,m7,m8,m9)

acf(resid(m9))
durbinWatsonTest(resid(m9))


pred <- jtools::make_predictions(m9,
                                 pred="EV",modx="Replay_differential",mod2="PCA1",
                                 at=list(Replay_differential=quantile(md$Replay_differential,c(.05,.95)),
                                         PCA1=quantile(md$PCA1,c(.05,.95))))
write.csv(pred,"ev_replay_PCA1.csv")

pred <- jtools::make_predictions(m9,
                                 pred="EV",modx="Replay_differential",mod2="PCA1",
                                 at=list(Replay_differential=quantile(md$Replay_differential,c(.05,.95)),
                                         PCA2=quantile(md$PCA2,c(.05,.95))))
write.csv(pred,"ev_replay_PCA2.csv")





interact_plot(m9,data=md,
              pred=EV,
              modx=Replay_differential,modx.values=quantile(d$Replay_differential,c(.05,.95)),
              mod2=PCA1,mod2.values=quantile(md$PCA1,c(.05,.95)),
              interval=TRUE,colors=c("#FF0055","#00FF8F"),vary.lty=FALSE) + 
  theme_classic()

interact_plot(m9,data=md,
              pred=EV,
              modx=Replay_differential,modx.values=quantile(d$Replay_differential,c(.05,.95)),
              mod2=PCA2,mod2.values=quantile(md$PCA2,c(.05,.95)),
              interval=TRUE,colors=c("#FF0055","#00FF8F"),vary.lty=FALSE) + 
  theme_classic()

g <- interact_plot(m9,data=md,
              pred=EV,
              mod2=Replay_differential,mod2.values=quantile(d$Replay_differential,c(.05,.95)),
              modx=PCA1,modx.values=unique(md$PCA1),
              interval=FALSE,vary.lty=FALSE) + 
  theme_classic() + scale_color_viridis(option="B",direction=-1)
ggsave("ev_replay_riskaversion.svg",g)

g <- interact_plot(m9,data=md,
              pred=EV,
              mod2=Replay_differential,mod2.values=quantile(d$Replay_differential,c(.05,.95)),
              modx=PCA2,modx.values=unique(md$PCA2),
              interval=FALSE,vary.lty=FALSE) + 
  theme_classic() + scale_color_viridis(option="B",direction=-1)
ggsave("ev_replay_anxiety.svg",g)



pd <- md %>%
  mutate(Choice=as.factor(Choice)) %>%
  group_by(Subject,Choice) %>%
  summarise(Replay=mean(Replay_differential),
            PCA1=mean(PCA1),PCA2=mean(PCA2)) %>%
  ungroup() %>%
  mutate(HighAnx = PCA2<median(PCA2),
         HighRiskAverse = PCA1<median(PCA1))


p1 <- ggplot(data=filter(pd,HighAnx==FALSE),aes(x=Choice,y=Replay,group=Subject,color=Choice)) + 
  geom_line(alpha=.5,colour="black") + 
  geom_point(size=3) + 
  theme_classic()
p2 <- ggplot(data=filter(pd,HighAnx==TRUE),aes(x=Choice,y=Replay,group=Subject,color=Choice)) + 
  geom_line(alpha=.5,colour="black") + 
  geom_point(size=3) + 
  theme_classic()
p3 <- ggplot(data=filter(pd,HighRiskAverse==FALSE),aes(x=Choice,y=Replay,group=Subject,color=Choice)) + 
  geom_line(alpha=.5,colour="black") + 
  geom_point(size=3) + 
  theme_classic()
p4 <- ggplot(data=filter(pd,HighRiskAverse==TRUE),aes(x=Choice,y=Replay,group=Subject,color=Choice)) + 
  geom_line(alpha=.5,colour="black") + 
  geom_point(size=3) + 
  theme_classic()

grid.arrange(p1,p2,p3,p4,nrow=2)

write.csv(pd,"choice_replay_personality.csv")




pred <- jtools::make_predictions(m9,
                                 pred="EV",modx="Replay_differential",
                                 at=list(
                                   Replay_differential=quantile(d$Replay_differential,c(.05,.95)),
                                   PCA1=quantile(md$PCA1,.05)))
pred <- evplots(pred)
write.csv(pred,"ev_replay_risk-high.csv")

pred <- jtools::make_predictions(m9,
                                 pred="EV",modx="Replay_differential",
                                 at=list(
                                   Replay_differential=quantile(d$Replay_differential,c(.05,.95)),
                                   PCA1=quantile(md$PCA1,.95)))
pred <- evplots(pred)
write.csv(pred,"ev_replay_risk-low.csv")

pred <- jtools::make_predictions(m9,
                                 pred="EV",modx="Replay_differential",
                                 at=list(
                                   Replay_differential=quantile(d$Replay_differential,c(.05,.95)),
                                   PCA2=quantile(md$PCA2,.05)))
pred <- evplots(pred)
write.csv(pred,"ev_replay_anx-high.csv")

pred <- jtools::make_predictions(m9,
                                 pred="EV",modx="Replay_differential",
                                 at=list(
                                   Replay_differential=quantile(d$Replay_differential,c(.05,.95)),
                                   PCA2=quantile(md$PCA2,.95)))
pred <- evplots(pred)
write.csv(pred,"ev_replay_anx-low.csv")


#######################################
# TEMPORAL MODULATION
#######################################

tm <- read.csv("D:/2020_RiskyReplay/results/replay/replay_temporal_lme.csv")
tm <- tm %>%
  filter(Subject!=680913,Subject!=263098)
tm <- cleandata(tm)
tm <- addQ(tm,Q)
tm$Modulation_average <- rowMeans(select(tm,Modulation_rewarding,Modulation_aversive))

# Add the two together
combined <- d %>% 
  filter(Lag>min(tm$Lag)-10,Lag<max(tm$Lag)+10) %>%
  mutate(Modulation_rewarding=tm$Modulation_rewarding,
         Modulation_aversive=tm$Modulation_aversive,
         Modulation_average=tm$Modulation_average) %>%
  mutate(Modulation_differential=Modulation_rewarding-Modulation_aversive)

md <- combined
md$EV <- scale(md$EV,center=TRUE,scale=FALSE)
md$Replay_differential <- scale(md$Replay_differential,center=TRUE,scale=FALSE)
md$Modulation_differential <- scale(md$Modulation_differential,center=TRUE,scale=FALSE)
md$Modulation_rewarding <- scale(md$Modulation_rewarding,center=TRUE,scale=FALSE)
md$Modulation_aversive <- scale(md$Modulation_aversive,center=TRUE,scale=FALSE)
md$BinModulation_differential <- md$Modulation_differential > 0
md$BinModulation_aversive <- md$Modulation_aversive > 0
md$BinModulation_rewarding <- md$Modulation_rewarding > 0

m0 <- glmer(Choice ~ EV*Replay_differential*Modulation_differential + RT + Certainty + (1|Subject/Lag),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

summary(m0)
interact_plot(m0,data=md,
              pred=EV,
              modx=Replay_differential,modx.values=quantile(d$Replay_differential,c(.05,.95)),
              mod2=Modulation_differential,mod2.values=quantile(md$Modulation_differential,c(.05,.95)),
              interval=TRUE,colors=c("#FF0055","#00FF8F"),vary.lty=FALSE) + 
  theme_classic()

pred <- jtools::make_predictions(m0,data=md,
              pred="EV",
              modx="Replay_differential",mod2="Modulation_differential",
              at=list(Replay_differential=quantile(d$Replay_differential,c(.05,.95)),
                      Modulation_differential=quantile(md$Modulation_differential,c(.05,.95)))) %>%
  group_by(EV) %>%
      mutate(grp=1:n())
pred$Replay_differential[pred$Replay_differential>0] = "rewarding"
pred$Replay_differential[pred$Replay_differential!="rewarding"] = "aversive"
pred$Modulation_differential[pred$Modulation_differential>0] = "rewarding_increase"
pred$Modulation_differential[pred$Modulation_differential!="rewarding_increase"] = "aversive_increase"

ggplot(data=pred,aes(x=EV,y=Choice,group=grp,color=Replay_differential)) + 
  geom_line(aes(linetype=Modulation_differential),size=1) + 
  theme_classic()
  

m1 <- glmer(Choice ~ EV*Replay_differential + EV*Modulation_rewarding + RT + Certainty + (1|Subject/Lag),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

summary(m1)
interact_plot(m1,data=md,
              pred=EV,
              modx=Replay_differential,modx.values=quantile(md$Replay_differential,c(.05,.95)),
              interval=TRUE,colors=c("#FF0055","#00FF8F"),vary.lty=FALSE) + 
  theme_classic()
interact_plot(m1,data=md,
              pred=EV,
              modx=Modulation_differential,modx.values=quantile(md$Modulation_differential,c(.05,.95)),
              interval=TRUE,colors=c("#FF0055","#00FF8F"),vary.lty=FALSE) + 
  theme_classic()


#######################################
# POST-PLANNING REPLAY
#######################################

postd <- read.csv("D:/2020_RiskyReplay/results/replay/replay_postplanning_lme.csv")
postd <- postd %>%
  filter(Subject!=680913,Subject!=263098)

postd <- cleandata(postd)
postd <- addQ(postd,Q)

postd$OutcomeType <- rep(NA,nrow(postd))
postd$OutcomeType[postd$Outcome<0] <- "bad"
postd$OutcomeType[postd$Outcome>0] <- "good"
postd$OutcomeType[postd$Choice==0] <- "safe"
postd$OutcomeType <- as.factor(postd$OutcomeType)

postd$ChoiceType <- rep(NA,nrow(postd))
postd$ChoiceType[postd$EV>0] <- "good"
postd$ChoiceType[postd$EV<0] <- "bad"
postd$ChoiceType <- as.factor(postd$ChoiceType)

postd$PE <- postd$Outcome - postd$EV

# add differential replay
thisd <- filter(d,Lag>min(postd$Lag)-10,Lag<max(postd$Lag)+10)
postd$Planning_differential <- rep(NA,nrow(postd))
for (subject in unique(postd$Subject)) {
  for (exptrial in unique(postd$ExpTrial[postd$Subject==subject])) {
    postd$Planning_differential[postd$Subject==subject & postd$ExpTrial==exptrial] <- thisd$Replay_differential[thisd$Subject==subject & thisd$ExpTrial==exptrial]
  }
}


md <- postd %>%
  filter(Lag>40,Lag<80)

md$Replay_differential <- scale(md$Replay_differential,center=TRUE,scale=FALSE)
md$Planning_differential <- scale(md$Planning_differential,center=TRUE,scale=FALSE)
md$Replay_average <- scale(md$Replay_average,center=TRUE,scale=FALSE)
md$Outcome <- scale(md$Outcome,center=TRUE,scale=FALSE)
md$EV <- scale(md$EV,center=TRUE,scale=FALSE)
md$PE <- scale(md$PE,center=TRUE,scale=FALSE)
md$Choice <- as.factor(md$Choice)

m0 <- lmer(Replay_differential ~ Outcome + RT + (1|Subject/Lag),
            data=filter(md,Choice==1), REML=TRUE,
            control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m1 <- lmer(Replay_differential ~ PE + RT + (1|Subject/Lag),
           data=filter(md,Choice==1), REML=TRUE,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m2 <- lmer(Replay_differential ~ EV*Outcome + RT + (1|Subject/Lag),
           data=filter(md,Choice==1), REML=TRUE,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m3 <- lmer(Replay_differential ~ Planning_differential + RT + (1|Subject/Lag),
           data=filter(md,Choice==1), REML=TRUE,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m4 <- lmer(Replay_differential ~ Outcome*Planning_differential + RT + (1|Subject/Lag),
           data=filter(md,Choice==1), REML=TRUE,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m5 <- lmer(Replay_differential ~ PE*Planning_differential + RT + (1|Subject/Lag),
           data=filter(md,Choice==1), REML=TRUE,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m6 <- lmer(Replay_differential ~ EV*Outcome*Planning_differential + RT + (1|Subject/Lag),
           data=filter(md,Choice==1), REML=TRUE,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

interact_plot(m5,pred=EV,modx=Outcome,modx.values=quantile(md$Outcome,c(.05,.95)),
              mod2=Planning_differential,mod2.values=quantile(md$Planning_differential,c(.05,.95)),
              interval=TRUE,colors=c("#FF0055","#00FF8F"),vary.lty=FALSE)


