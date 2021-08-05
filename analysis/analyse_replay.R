library(lme4)
library(lmerTest)
library(interactions)
library(car)
library(emmeans)
library(tidyverse)
library(viridis)
library(extrafont)

#######################################
# LOAD DATA
#######################################

d <- read.csv("D:/2020_RiskyReplay/results/replay/replay_lme.csv")
Q <- read.csv("D:/2020_RiskyReplay/results/questionnaire_results.csv")

#######################################
# CLEAN DATA
#######################################

d$Subject <- as.factor(d$Subject)
# d$Choice <- as.factor(d$Choice) # 1 is approach, 0 is avoid
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

d$Replay_average <- (d$Replay_rewarding + d$Replay_aversive)/2

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
library(ggcorrplot)
C <- cbind(Q,PC) %>% select(starts_with("score_"),starts_with("PC",ignore.case=FALSE))
corr <- round(cor(C),1)
p.mat <- cor_pmat(C)
ggcorrplot(corr,type="upper",lab=FALSE)


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

d$IUS <- scale(d$IUS,center=TRUE,scale=FALSE)
d$Worry <- scale(d$Worry,center=TRUE,scale=FALSE)
d$Risk <- scale(d$Risk,center=TRUE,scale=FALSE)
d$PCA1 <- scale(d$PCA1,center=TRUE,scale=FALSE)
d$PCA2 <- scale(d$PCA2,center=TRUE,scale=FALSE)

d$Anxiety <- rowMeans(select(d,IUS,Worry))

d$OverallBlock <- rep(NA,nrow(d))
d$OverallBlock[d$Practice==1] <- 0
for (b in 1:10) {
  d$OverallBlock[d$Practice==0 & d$Block==b] <- b
}

#######################################
# BEHAVIOURAL ANALYSIS
#######################################

md <- filter(d,Lag==d$Lag[1])
# md$EV <- scale(md$EV,center=TRUE,scale=FALSE)
md$PCA1 <- scale(md$PCA1,center=TRUE,scale=FALSE)
md$PCA2 <- scale(md$PCA2,center=TRUE,scale=FALSE)

m0 <- glmer(Choice ~ EV + (1|Subject),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m1 <- glmer(Choice ~ EV + Certainty + (1|Subject),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m2 <- glmer(Choice ~ EV*Certainty + (1|Subject),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

anova(m0,m1,m2)

m3 <- glmer(Choice ~ EV*PCA1 + Certainty + (1|Subject),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m4 <- glmer(Choice ~ EV*PCA2 + Certainty + (1|Subject),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m5 <- glmer(Choice ~ EV*PCA1 + EV*PCA2 + Certainty + (1|Subject),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

anova(m2,m3,m4,m5)

pred <- jtools::make_predictions(m0,pred="EV")

ggplot(data=pred, aes(x=EV,y=Choice)) + 
  geom_ribbon(aes(ymax=ymax,ymin=ymin),alpha=.2) + 
  geom_line(size=1) +
  theme_classic() +
  coord_cartesian(xlim=c(-9,11)) + scale_x_continuous(breaks=seq(-9,11,2))

#######################################
# REPLAY ANALYSIS (VALUE)
#######################################

md <- d #filter(d,Subject!=97403)

md <- md %>%
  group_by(Subject,OverallBlock) %>%
  mutate(BlockAcc=mean(Acc=="correct")) %>%
  group_by(Subject) %>%
  mutate(SubjAcc=mean(Acc=="correct"))

maxtype <- "group" # subject or group

if (maxtype=="subject"){
  
  # Choose subject's absolute maximum
  bestlags <- md %>% 
    filter(Lag>30 & Lag<90) %>%
    group_by(Subject,Lag) %>%
    summarise(Replay_average=mean(Replay_average)) %>%
    group_by(Subject) %>%
    slice(which.max(abs(Replay_average))) %>%
    as.data.frame()
  
  for (subject in unique(md$Subject)) {
    
    thislag <- bestlags$Lag[bestlags$Subject==subject]
    thisd <- filter(md,Subject==subject,Lag==thislag)
    thisd$BestLag <- rep(thislag,nrow(thisd))
    if (subject==unique(md$Subject)[1]){
      newmd <- thisd
    } else {
      newmd <- rbind(newmd,thisd)
    }
    
  }
  md <- newmd
  
} else if (maxtype=="group"){
  
  # # OR choose the overall maximum (60 ms)
  # md <- md %>% filter(Lag==60)
  
  # OR average across a range of lags
  md <- md %>% 
    filter(Lag>10 & Lag<100) %>%
    group_by(Subject,ExpTrial) %>%
    mutate(Replay_differential=mean(Replay_differential),
           Replay_rewarding=mean(Replay_rewarding),
           Replay_aversive=mean(Replay_aversive)) %>%
    filter(Lag==60)
  
}

# Rearrange so that replay is spread out across rows instead of columns
md <- md %>%
  gather(Replay_type, Sequenceness, c(Replay_rewarding:Replay_differential,Replay_average), factor_key=TRUE) %>%
  mutate(Replay_type=gsub("Replay_","",Replay_type))

md$PathVal <- rep(NA,nrow(md))
md$PathVal[md$Replay_type=="rewarding"] <- apply(select(filter(md,Replay_type=="rewarding"),nV_1:nV_1),1,max)
md$PathVal[md$Replay_type=="aversive"] <- apply(select(filter(md,Replay_type=="aversive"),nV_1:nV_1),1,min)

md$PathProb <- rep(NA,nrow(md))
md$PathProb[md$Replay_type=="rewarding"] <- select(filter(md,Replay_type=="rewarding"),RewProb)
md$PathProb[md$Replay_type=="aversive"] <- 1-select(filter(md,Replay_type=="rewarding"),RewProb)

# Make models

m0 <- lmer(Sequenceness ~ PathVal + (1|Subject),
            data=filter(md,Replay_type=="rewarding" | Replay_type=="aversive"), REML=TRUE,
            control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))


#######################################
# REPLAY ANALYSIS (CHOICE)
#######################################

md <- d #filter(d,Subject!=97403)

md <- md %>%
  group_by(Subject,OverallBlock) %>%
  mutate(BlockAcc=mean(Acc=="correct")) %>%
  group_by(Subject) %>%
  mutate(SubjAcc=mean(Acc=="correct"))

maxtype <- "group" # subject or group

if (maxtype=="subject"){
  
  # Choose subject's absolute maximum
  bestlags <- md %>% 
    filter(Lag>30 & Lag<90) %>%
    group_by(Subject,Lag) %>%
    summarise(Replay_average=mean(Replay_average)) %>%
    group_by(Subject) %>%
    slice(which.max(abs(Replay_average))) %>%
    as.data.frame()
  
  for (subject in unique(md$Subject)) {
    
    thislag <- bestlags$Lag[bestlags$Subject==subject]
    thisd <- filter(md,Subject==subject,Lag==thislag)
    thisd$BestLag <- rep(thislag,nrow(thisd))
    if (subject==unique(md$Subject)[1]){
      newmd <- thisd
    } else {
      newmd <- rbind(newmd,thisd)
    }
    
  }
  md <- newmd
  
} else if (maxtype=="group"){
  
  # # OR choose the overall maximum (60 ms)
  # md <- md %>% filter(Lag==60)
  
  # OR average across a range of lags
  md <- md %>% 
    filter(Lag>10 & Lag<100) %>%
    group_by(Subject,ExpTrial) %>%
    mutate(Replay_differential=mean(Replay_differential)) %>%
    filter(Lag==60)
  
}

# # remove outliers
# md <- md %>% 
#   as.data.frame() %>%
#   # group_by(Subject) %>%
#   mutate(outliers=abs(scale(Replay_differential))>3) %>%
#   filter(outliers==FALSE)

# scale_type <- "group" # whether to scale within subject or just mean-centre the group
# if (scale_type=="subject"){
#   md <- md %>%
#     group_by(Subject) %>%
#     mutate(EV=scale(EV),Replay_differential=scale(Replay_differential))
# } else if (scale_type=="group") {
  md$EV <- scale(md$EV,center=TRUE,scale=FALSE)
  md$Replay_differential <- scale(md$Replay_differential,center=TRUE,scale=FALSE)
# }

m0 <- glmer(Choice ~ EV + (1|Subject),
           data=md, family="binomial",
           control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m1 <- glmer(Choice ~ EV*Replay_differential + (1|Subject),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m2 <- glmer(Choice ~ EV*Replay_differential*PCA1 + (1|Subject),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

m3 <- glmer(Choice ~ EV*Replay_differential*PCA2 + (1|Subject),
           data=md, family="binomial",
           control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

anova(m0,m1,m2,m3)




interact_plot(m3,data=md,pred=EV,modx=Replay_differential,modx.values=quantile(d$Replay_differential,c(.05,.5,.95)),
              interval=TRUE,colors=c("#FF0055","#737373","#00FF8F"),vary.lty=FALSE) + 
  theme_classic()
