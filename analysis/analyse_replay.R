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
library(ggcorrplot)

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

# scree plot
screeplot(pca,type="l")
abline(h=1,col="red",lty=5)

# correlation matrix
C <- cbind(Q,PC) %>% select(starts_with("score_"),starts_with("PC",ignore.case=FALSE))
corr <- round(cor(C),1)
p.mat <- cor_pmat(C)
ggcorrplot(corr,type="upper",lab=FALSE)


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

# Plot correlations
pd <- cbind(Q,PC) %>%
  gather(PCAnum,PCAval,PC1:PC2)

p1 <- ggplot(data=pd,aes(x=score_IUS,y=PCAval,group=PCAnum)) +
  geom_point() +
  geom_smooth(method=lm) +
  facet_wrap(~PCAnum) + theme_classic() +
  coord_cartesian(xlim=c(10,55)) + scale_x_continuous(breaks=seq(10,55,5))

ggplot(data=pd,aes(x=score_worry,y=PCAval,group=PCAnum)) +
  geom_point() +
  geom_smooth(method=lm) +
  facet_wrap(~PCAnum) + theme_classic()

ggplot(data=pd,aes(x=score_risk,y=PCAval,group=PCAnum)) +
  geom_point() +
  geom_smooth(method=lm) +
  facet_wrap(~PCAnum) + theme_classic()

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

# MODEL 1
m1 <- glmer(Choice ~ RewVal*LossVal*RewProb + Certainty + RT + (1|Subject),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

summary(m1)
round(range(vif(m1)),2)
round(durbinWatsonTest(resid(m1)),2)




p1 <- interact_plot(m1,pred=RewProb,modx=RewVal,modx.values=quantile(md$RewVal,c(.05,.5,.95)),
              interval=TRUE,vary.lty=FALSE) + 
  theme_classic()

p2 <- interact_plot(m1,pred=RewProb,modx=LossVal,modx.values=quantile(md$LossVal,c(.05,.5,.95)),
                    interval=TRUE,vary.lty=FALSE) + 
  theme_classic()

g <- grid.arrange(p1,p2,nrow=1)

ggsave("rewvalrewprob.svg",g)


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

md$PathRecency <- rep(NA,nrow(md))
md$PathRecency[md$Replay_type=="rewarding"] <- md$RewPathRecency_block[md$Replay_type=="rewarding"]
md$PathRecency[md$Replay_type=="aversive"] <- md$LossPathRecency_block[md$Replay_type=="aversive"]
md$PathRecency <- log(md$PathRecency+1)

md$Choice <- as.factor(md$Choice)

md <- filter(md,
             Lag>10,Lag<100,
             Replay_type=="rewarding" | Replay_type=="aversive")

md$PathExp <- scale(md$PathExp,center=TRUE,scale=FALSE)
md$PathVal <- scale(md$PathVal,center=TRUE,scale=FALSE)
md$RT <- scale(md$RT,center=TRUE,scale=FALSE)
md$PathProb <- scale(md$PathProb,center=TRUE,scale=FALSE)
md$PathRecency <- scale(md$PathRecency,center=TRUE,scale=FALSE)

# MODEL 2
m2 <- lmer(Sequenceness ~ Replay_type*Choice + RT + (1|Subject/Lag),
           data=md, REML=FALSE,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

summary(m2)
round(range(vif(m2)),2)
round(durbinWatsonTest(resid(m2)),2)


# MODEL 3
m3 <- lmer(Sequenceness ~ PathRecency*Replay_type*PathProb + RT + (1|Subject/Lag),
           data=md, REML=TRUE,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))
summary(m3)
round(range(vif(m3)),2)
round(durbinWatsonTest(resid(m3)),2)


pred <- jtools::make_predictions(m3,pred="PathRecency",
                                 at=list(Replay_type=unique(md$Replay_type),
                                         PathProb=unique(md$PathProb))) %>%
  mutate(semupper=Sequenceness+((ymax-ymin)/3.96),
         semlower=Sequenceness-((ymax-ymin)/3.96))

ggplot(pred,aes(x=PathRecency,y=Sequenceness,
                group=Replay_type,color=PathProb)) + 
  geom_point() + 
  theme_classic() + 
  coord_cartesian(ylim=c(0,0.04))


# MODEL 4
m4 <- lmer(Sequenceness ~ PathProb*Choice + RT + (1|Subject/Lag),
           data=md, REML=FALSE,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))
summary(m4)
round(range(vif(m4)),2)
round(durbinWatsonTest(resid(m4)),2)

pred <- jtools::make_predictions(m4,pred="PathProb",
                                 at=list(Choice=unique(md$Choice))) %>%
  mutate(semupper=Sequenceness+((ymax-ymin)/3.96),
         semlower=Sequenceness-((ymax-ymin)/3.96))

ggplot(pred,aes(x=PathProb,y=Sequenceness,group=Choice,color=Choice,fill=Choice)) + 
  geom_ribbon(aes(ymin=semlower,ymax=semupper),alpha=.2,color=NA) + 
  geom_line(size=1) + 
  theme_classic() + 
  coord_cartesian(ylim=c(-0,0.04))

#######################################
# REPLAY ANALYSIS (CHOICE)
#######################################

md <- d %>%
  filter(Lag>10, Lag<100) # Subject!=97403

md$origEV <- md$EV
md$EV <- scale(md$EV,center=TRUE,scale=FALSE)
md$RT <- scale(md$RT,center=TRUE,scale=FALSE)
md$Replay_differential <- scale(md$Replay_differential,center=TRUE,scale=FALSE)
md$Replay_rewarding <- scale(md$Replay_rewarding,center=TRUE,scale=FALSE)
md$Replay_aversive <- scale(md$Replay_aversive,center=TRUE,scale=FALSE)

# MODEL 5
m5 <- glmer(Choice ~ EV*Replay_differential + Certainty + RT + (1|Subject/Lag),
           data=md, family="binomial",
           control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))
summary(m5)
round(range(vif(m5)),2)
round(durbinWatsonTest(resid(m5)),2)

interact_plot(m5,data=md,pred=EV,modx=Replay_differential,modx.values=quantile(d$Replay_differential,c(.05,.5,.95)),
              interval=TRUE,colors=c("#FF0055","#737373","#00FF8F"),vary.lty=FALSE) + 
  theme_classic()

# MODEL 6
m6 <- glmer(Choice ~ EV*Replay_rewarding+ EV*Replay_aversive + Certainty + RT + (1|Subject/Lag),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))
summary(m6)
round(range(vif(m6)),2)
round(durbinWatsonTest(resid(m6)),2)

interact_plot(m6,data=md,pred=EV,modx=Replay_rewarding,modx.values=quantile(d$Replay_rewarding,c(.05,.95)),
              interval=TRUE,vary.lty=FALSE) + 
  theme_classic()
interact_plot(m6,data=md,pred=EV,modx=Replay_aversive,modx.values=quantile(d$Replay_aversive,c(.05,.95)),
              interval=TRUE,vary.lty=FALSE) + 
  theme_classic()


# Data to make plots for model 5 and model 6
pred <- jtools::make_predictions(m5,
                                 pred="EV",modx="Replay_differential",
                                 at=list(Replay_differential=quantile(d$Replay_differential,c(.05,.5,.95))))
names(pred)[names(pred)=="Replay_differential"] <- "Replay"
write.csv(pred,"ev_replay_differential.csv")

pred <- jtools::make_predictions(m6,
                                 pred="EV",modx="Replay_rewarding",
                                 at=list(Replay_rewarding=quantile(d$Replay_rewarding,c(.05,.5,.95))))
names(pred)[names(pred)=="Replay_rewarding"] <- "Replay"
write.csv(pred,"ev_replay_reward.csv")

pred <- jtools::make_predictions(m6,
                                 pred="EV",modx="Replay_aversive",
                                 at=list(Replay_aversive=quantile(d$Replay_aversive,c(.05,.5,.95))))
names(pred)[names(pred)=="Replay_aversive"] <- "Replay"
write.csv(pred,"ev_replay_loss.csv")


ss <- sim_slopes(m5,pred=EV,modx=Replay_differential,modx.values=quantile(d$Replay_differential,c(.05,.95)),
                 cond.int=TRUE)

# MODEL 7
m7 <- glmer(Choice ~ EV*Replay_differential*PCA1 + EV*Replay_differential*PCA2 + Certainty + RT + (1|Subject/Lag),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))
summary(m7)
round(range(vif(m7)),2)
round(durbinWatsonTest(resid(m7)),2)

# MODEL 8
m8 <- glmer(Choice ~ EV*Replay_rewarding*PCA1 + EV*Replay_rewarding*PCA2 + 
              EV*Replay_aversive*PCA1 + EV*Replay_aversive*PCA2 + 
              Certainty + RT + (1|Subject/Lag),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))
summary(m8)
round(range(vif(m8)),2)
round(durbinWatsonTest(resid(m8)),2)

# MODEL COMPARISONS FOR ANXIETY AND RISK AVERSION
mA <- glmer(Choice ~ EV*Replay_differential + Certainty + RT + (1|Subject/Lag),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

mB <- glmer(Choice ~ EV*Replay_differential + PCA1 + Certainty + (1|Subject/Lag),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

mC <- glmer(Choice ~ EV*Replay_differential + PCA2 + Certainty + (1|Subject/Lag),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

mD <- glmer(Choice ~ EV*Replay_differential + PCA1 + PCA2 + Certainty + (1|Subject/Lag),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

mE <- glmer(Choice ~ EV*Replay_differential*PCA1 + Certainty + (1|Subject/Lag),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

mF <- glmer(Choice ~ EV*Replay_differential*PCA1 + PCA2 + Certainty + (1|Subject/Lag),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

mG <- glmer(Choice ~ EV*Replay_differential*PCA2 + Certainty + (1|Subject/Lag),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

mH <- glmer(Choice ~ EV*Replay_differential*PCA2 + PCA1 + Certainty + (1|Subject/Lag),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

mI <- glmer(Choice ~ EV*Replay_differential*PCA1 + EV*Replay_differential*PCA2 + Certainty + (1|Subject/Lag),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

anova(mA,mB,mC,mD,mE,mF,mG,mH,mI)


# plots for model 7
interact_plot(m7,data=md,
              pred=EV,
              modx=Replay_differential,modx.values=quantile(d$Replay_differential,c(.05,.95)),
              mod2=PCA1,mod2.values=quantile(md$PCA1,c(.05,.95)),
              interval=TRUE,colors=c("#FF0055","#00FF8F"),vary.lty=FALSE) + 
  theme_classic()

interact_plot(m7,data=md,
              pred=EV,
              modx=Replay_differential,modx.values=quantile(d$Replay_differential,c(.05,.95)),
              mod2=PCA2,mod2.values=quantile(md$PCA2,c(.05,.95)),
              interval=TRUE,colors=c("#FF0055","#00FF8F"),vary.lty=FALSE) + 
  theme_classic()

g <- interact_plot(m7,data=md,
              pred=EV,
              mod2=Replay_differential,mod2.values=quantile(d$Replay_differential,c(.05,.95)),
              modx=PCA1,modx.values=unique(md$PCA1),
              interval=FALSE,vary.lty=FALSE) + 
  theme_classic() + scale_color_viridis(option="B",direction=-1)
ggsave("ev_replay_riskaversion.svg",g)

g <- interact_plot(m7,data=md,
              pred=EV,
              mod2=Replay_differential,mod2.values=quantile(d$Replay_differential,c(.05,.95)),
              modx=PCA2,modx.values=unique(md$PCA2),
              interval=FALSE,vary.lty=FALSE) + 
  theme_classic() + scale_color_viridis(option="B",direction=-1)
ggsave("ev_replay_anxiety.svg",g)




pred <- jtools::make_predictions(m7,
                                 pred="EV",modx="Replay_differential",mod2="PCA1",
                                 at=list(Replay_differential=quantile(md$Replay_differential,c(.05,.95)),
                                         PCA1=quantile(md$PCA1,c(.05,.95))))
write.csv(pred,"ev_replay_PCA1.csv")

pred <- jtools::make_predictions(m7,
                                 pred="EV",modx="Replay_differential",mod2="PCA1",
                                 at=list(Replay_differential=quantile(md$Replay_differential,c(.05,.95)),
                                         PCA2=quantile(md$PCA2,c(.05,.95))))
write.csv(pred,"ev_replay_PCA2.csv")

pred <- jtools::make_predictions(m7,
                                 pred="EV",modx="Replay_differential",
                                 at=list(
                                   Replay_differential=quantile(d$Replay_differential,c(.05,.95)),
                                   PCA1=quantile(md$PCA1,.05)))
pred <- evplots(pred)
write.csv(pred,"ev_replay_risk-high.csv")

pred <- jtools::make_predictions(m7,
                                 pred="EV",modx="Replay_differential",
                                 at=list(
                                   Replay_differential=quantile(d$Replay_differential,c(.05,.95)),
                                   PCA1=quantile(md$PCA1,.95)))
pred <- evplots(pred)
write.csv(pred,"ev_replay_risk-low.csv")

pred <- jtools::make_predictions(m7,
                                 pred="EV",modx="Replay_differential",
                                 at=list(
                                   Replay_differential=quantile(d$Replay_differential,c(.05,.95)),
                                   PCA2=quantile(md$PCA2,.05)))
pred <- evplots(pred)
write.csv(pred,"ev_replay_anx-high.csv")

pred <- jtools::make_predictions(m7,
                                 pred="EV",modx="Replay_differential",
                                 at=list(
                                   Replay_differential=quantile(d$Replay_differential,c(.05,.95)),
                                   PCA2=quantile(md$PCA2,.95)))
pred <- evplots(pred)
write.csv(pred,"ev_replay_anx-low.csv")

