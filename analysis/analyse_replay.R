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
R <- read.csv("D:/2020_RiskyReplay/results/reactivation/reactivation.csv")

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
# ADD THE OVERALL REACTIVATION
#######################################

d$Ymult_rewarding <- rep(NA,nrow(d))
d$Ymult_aversive <- rep(NA,nrow(d))
d$Ysig_rewarding <- rep(NA,nrow(d))
d$Ysig_aversive <- rep(NA,nrow(d))

for (subject in unique(d$Subject)){
  
  print(subject)

  for (trial in unique(d$ExpTrial[d$Subject==subject])){
    
    idx <- d$Subject==subject & d$ExpTrial==trial
    r <- R[R$Subject==subject & R$ExpTrial==trial,]
    
    d$Ymult_rewarding[idx] <- r$Ymult_reactivation[r$Reactivated_path=="rewarding"]
    d$Ymult_aversive[idx] <- r$Ymult_reactivation[r$Reactivated_path=="aversive"]
    d$Ysig_rewarding[idx] <- r$Ysig_reactivation[r$Reactivated_path=="rewarding"]
    d$Ysig_aversive[idx] <- r$Ysig_reactivation[r$Reactivated_path=="aversive"]
    
    d$Mmult_rewarding[idx] <- r$Mmult_reactivation[r$Reactivated_path=="rewarding"]
    d$Mmult_aversive[idx] <- r$Mmult_reactivation[r$Reactivated_path=="aversive"]
    d$Msig_rewarding[idx] <- r$Msig_reactivation[r$Reactivated_path=="rewarding"]
    d$Msig_aversive[idx] <- r$Msig_reactivation[r$Reactivated_path=="aversive"]
    
  }
}

#########################################################
# ADD BEHAVIOURAL STRATEGIES & COUNTERFACTUAL UTILITY
#########################################################

bd <- read.csv('D:/2020_RiskyReplay/results/modelling/behavmodeltable.csv')
countd <- read.csv('D:/2020_RiskyReplay/results/modelling/counterfactualtable.csv')

d$bModelNum <- rep(NA,nrow(d))
d$bCalculationType <- rep(NA,nrow(d))
d$bPathChoice <- rep(NA,nrow(d))

d$utility_goal_mCalculate <- rep(NA,nrow(d))
d$utility_counter_mCalculate <- rep(NA,nrow(d))
d$utility_goal_mLearn <- rep(NA,nrow(d))
d$utility_counter_mLearn <- rep(NA,nrow(d))

for (subject in unique(d$Subject)){
  
  d$bModelNum[d$Subject==subject] = bd$bModelNum[bd$Subject==subject]
  d$bCalculationType[d$Subject==subject] = bd$bCalculationType[bd$Subject==subject]
  d$bPathChoice[d$Subject==subject] = bd$bPathChoice[bd$Subject==subject]
  
  d$utility_goal_mCalculate[d$Subject==subject] = countd$utility_goal_mCalculate[countd$Subject==subject]
  d$utility_counter_mCalculate[d$Subject==subject] = countd$utility_counter_mCalculate[countd$Subject==subject]
  d$utility_goal_mLearn[d$Subject==subject] = countd$utility_goal_mLearn[countd$Subject==subject]
  d$utility_counter_mLearn[d$Subject==subject] = countd$utility_counter_mLearn[countd$Subject==subject]
  
}

d$bModelNum <- as.factor(d$bModelNum)
d$bCalculationType <- as.factor(d$bCalculationType)
d$bPathChoice <- as.factor(d$bPathChoice)


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
# POLICY CHANGE
#######################################

d$pchange_EV <- rep(NA,nrow(d)) # policy change at next same EV
d$pchange_EV_distance <- rep(NA,nrow(d))
d$pchange_next <- rep(NA,nrow(d))
d$pchange_nextEV <- rep(NA,nrow(d))
d$pchange_nextAcc <- rep(NA,nrow(d))
d$pchange_thisImprovement <- rep(NA,nrow(d))
d$pchange_nextImprovement <- rep(NA,nrow(d))
for (subject in unique(d$Subject)){
  
  thisd <- filter(d,Subject==subject)
  
  trials <- unique(thisd$ExpTrial)
  
  for (trl in 1:length(trials)){
    
    currentChoice <- thisd$Choice[thisd$ExpTrial==trials[trl]][1] # take the first, as there will be multiple lags
    nextChoice <- thisd$Choice[thisd$ExpTrial==trials[trl+1]][1]
    
    currentEV <- round(thisd$EV[thisd$ExpTrial==trials[trl]][1])
    nextEV <- which(round(thisd$EV)==currentEV & thisd$ExpTrial > thisd$ExpTrial[trl])
    if (length(nextEV)>0){
      nextEVchoice <- thisd$Choice[thisd$ExpTrial[nextEV[1]]][1]
      nextEVdist <- thisd$ExpTrial[nextEV[1]] - trials[trl]
    } else {
      nextEVchoice <- NA
      nextEVdist <- NA
    }
    
    idx <- d$Subject==subject & d$ExpTrial==trials[trl]
    
    d$pchange_EV[idx] <- currentChoice != nextEVchoice
    d$pchange_EV_distance[idx] <- nextEVdist
    d$pchange_EV_nextAcc[idx] <- thisd$Acc[thisd$ExpTrial[nextEV[1]]][1]
    
    d$pchange_next[idx] <- currentChoice != nextChoice
    d$pchange_nextEV[idx] <- round(thisd$EV[thisd$ExpTrial==trials[trl+1]][1])
    
    d$pchange_nextAcc[idx] <- thisd$Acc[thisd$ExpTrial==trials[trl+1]][1]
    
    thisPerformance <- thisd$EV[thisd$ExpTrial==trials[trl]][1]
    if (currentChoice==0){
      thisPerformance <- (-1)*thisPerformance
    }
    nextPerformance <- thisd$EV[thisd$ExpTrial==trials[trl+1]][1]
    if (currentChoice==0){
      nextPerformance <- (-1)*nextPerformance
    }
    
    d$pchange_thisImprovement[idx] <- thisPerformance
    d$pchange_nextImprovement[idx] <- nextPerformance
    
  }
}
d$pchange_EV_nextAcc[d$pchange_EV_nextAcc==2] <- 0

# Creat new 'counterfactual' replay variable - i.e. how much more the SUBOPTIMAL OUTCOME of a planned choice is replayed
d$Replay_counterfactual <- rep(NA,nrow(d))
d$Replay_counterfactual[d$Choice==1] <- d$Replay_aversive[d$Choice==1] - d$Replay_rewarding[d$Choice==1]
d$Replay_counterfactual[d$Choice==0] <- d$Replay_rewarding[d$Choice==0] - d$Replay_aversive[d$Choice==0]

# Shift replay to NEXT trial
d$Replay_next_differential <- rep(NA,nrow(d))
for (subject in unique(d$Subject)){
  
  thisd <- filter(d,Subject==subject)
  
  trials <- unique(thisd$ExpTrial)
  
  for (trl in 1:(length(trials)-1)){

    idx1 <- d$Subject==subject & d$ExpTrial==trials[trl]
    idx2 <- d$Subject==subject & d$ExpTrial==trials[trl+1]
    
    d$Replay_next_differential[idx1] <- d$Replay_differential[idx2]
    
  }
}

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

md$PathNum <- rep(NA,nrow(md))
md$PathNum[md$Replay_type=="rewarding" & md$Path1Type=="rewarding"] <- 1
md$PathNum[md$Replay_type=="rewarding" & md$Path1Type=="aversive"] <- 2
md$PathNum[md$Replay_type=="aversive" & md$Path1Type=="rewarding"] <- 2
md$PathNum[md$Replay_type=="aversive" & md$Path1Type=="aversive"] <- 1
md$PathNum <- as.factor(md$PathNum)

md$PathVal <- rep(NA,nrow(md))
md$PathVal[md$Replay_type=="rewarding"] <- apply(select(filter(md,Replay_type=="rewarding"),nV_1:nV_2),1,max)
md$PathVal[md$Replay_type=="aversive"] <- apply(select(filter(md,Replay_type=="aversive"),nV_1:nV_2),1,min)

md$PathProb <- rep(NA,nrow(md))
md$PathProb[md$Replay_type=="rewarding"] <- md$RewProb[md$Replay_type=="rewarding"]
md$PathProb[md$Replay_type=="aversive"] <- round(1-md$RewProb[md$Replay_type=="aversive"],1)

md$PathCertainty <- md$PathProb
md$PathCertainty[which(md$PathProb==0.1)] <- rep(0.9,sum(md$PathProb==0.1,na.rm=TRUE))
md$PathCertainty[which(md$PathProb==0.3)] <- rep(0.7,sum(md$PathProb==0.3,na.rm=TRUE))

md$PathRecency <- rep(NA,nrow(md))
md$PathRecency[md$Replay_type=="rewarding"] <- md$RewPathRecency_block[md$Replay_type=="rewarding"]
md$PathRecency[md$Replay_type=="aversive"] <- md$LossPathRecency_block[md$Replay_type=="aversive"]
md$PathRecency <- log(md$PathRecency+1)

md$Choice <- as.factor(md$Choice)

md <- filter(md,
             Lag>10,Lag<100,
             Replay_type=="rewarding" | Replay_type=="aversive")

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

e <- emmeans(m2,pairwise~Replay_type*Choice)

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
md$Replay_counterfactual <- scale(md$Replay_counterfactual,center=TRUE,scale=FALSE)

md$utility_goal_mCalculate <- scale(md$utility_goal_mCalculate,center=TRUE,scale=FALSE)
md$utility_counter_mCalculate <- scale(md$utility_counter_mCalculate,center=TRUE,scale=FALSE)
md$utility_goal_mLearn <- scale(md$utility_goal_mLearn,center=TRUE,scale=FALSE)
md$utility_counter_mLearn <- scale(md$utility_counter_mLearn,center=TRUE,scale=FALSE)

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


# MODEL ???? - POLICY CHANGE
md <- d %>%
  filter(Lag>10, Lag<100, !is.na(pchange_EV), !is.na(pchange_next)) # Subject!=97403

md$PE <- md$Outcome - md$EV
md$PE <- scale(md$PE,center=TRUE,scale=FALSE)

md$origEV <- md$EV
md$EV <- scale(md$EV,center=TRUE,scale=FALSE)
md$RT <- scale(md$RT,center=TRUE,scale=FALSE)
md$Replay_differential <- scale(md$Replay_differential,center=TRUE,scale=FALSE)
md$Replay_rewarding <- scale(md$Replay_rewarding,center=TRUE,scale=FALSE)
md$Replay_aversive <- scale(md$Replay_aversive,center=TRUE,scale=FALSE)
md$Replay_counterfactual <- scale(md$Replay_counterfactual,center=TRUE,scale=FALSE)
md$pchange_EV_distance_exp <- scale(exp(md$pchange_EV_distance),center=TRUE,scale=FALSE)
md$pchange_EV_distance <- scale(md$pchange_EV_distance,center=TRUE,scale=FALSE)
md$pchange_nextEV <- scale(md$pchange_nextEV,center=TRUE,scale=FALSE)
md$Outcome <- scale(md$Outcome,center=TRUE,scale=FALSE)
md$pchange_nextImprovement <- scale(md$pchange_nextImprovement,center=TRUE,scale=FALSE)
md$pchange_thisImprovement <- scale(md$pchange_thisImprovement,center=TRUE,scale=FALSE)

md$Counterfactual_path_recency <- rep(NA,nrow(md))
md$Counterfactual_path_recency[md$Choice==1] <- md$RewPathExp_half[md$Choice==1]
md$Counterfactual_path_recency[md$Choice==0] <- md$LossPathExp_half[md$Choice==0]

md$pchange_EV[md$pchangeEV=="yes"] <- 1
md$pchange_EV[md$pchangeEV=="no"] <- 0
md$pchange_EV <- as.factor(md$pchange_EV)

md$pchange_next[md$pchangeEV=="yes"] <- 1
md$pchange_next[md$pchangeEV=="no"] <- 0
md$pchange_next <- as.factor(md$pchange_EV)

md$pchange_EV_nextAcc <- as.factor(md$pchange_EV_nextAcc)


m_EV <- glmer(pchange_EV ~ Choice*Replay_counterfactual*EV + pchange_EV_distance_exp + RT + (1|Subject/Lag),
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))
summary(m_EV)
round(range(vif(m_EV)),2)
round(durbinWatsonTest(resid(m_EV)),2)
interact_plot(m_EV,pred=EV,modx=Replay_counterfactual,mod2=Choice,
              interval=TRUE,vary.lty=FALSE)

# more complete model
m_next <- glmer(pchange_next ~ Replay_counterfactual*pchange_nextAcc*Acc + (1|Subject/Lag),
              data=md, family="binomial",
              control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))
summary(m_next)
round(range(vif(m_next)),2)
round(durbinWatsonTest(resid(m_next)),2)
interact_plot(m_next,pred=Replay_counterfactual,modx=Choice,mod2=EV,
              interval=TRUE,vary.lty=FALSE)


anova(m_EV,m_next)





pd <- md %>%
  group_by(Subject,Choice,pchange_next) %>%
  summarise(replay=mean(Replay_counterfactual))

pdsum <- pd %>%
  group_by(Choice,pchange_next) %>%
  summarise(replay=mean(replay))

ggplot(pdsum,aes(x=Choice,y=replay,group=pchange_next)) + geom_bar()

#######################################
# OVERALL STATE REACTIVATION
#######################################

# Rearrange so that replay is spread out across rows instead of columns
md <- d %>%
  gather(Replay_type, Sequenceness, c(Replay_rewarding:Replay_differential,Replay_average), factor_key=TRUE) %>%
  mutate(Replay_type=gsub("Replay_","",Replay_type))
md$Replay_type <- as.factor(md$Replay_type)

md$Ymult <- rep(NA,nrow(md))
md$Ymult[md$Replay_type=="rewarding"] <- md$Ymult_rewarding[md$Replay_type=="rewarding"]
md$Ymult[md$Replay_type=="aversive"] <- md$Ymult_aversive[md$Replay_type=="aversive"]
md$Ymult[md$Replay_type=="differential"] <- md$Ymult_rewarding[md$Replay_type=="differential"] - md$Ymult_aversive[md$Replay_type=="differential"]

md$Ysig <- rep(NA,nrow(md))
md$Ysig[md$Replay_type=="rewarding"] <- md$Ysig_rewarding[md$Replay_type=="rewarding"]
md$Ysig[md$Replay_type=="aversive"] <- md$Ysig_aversive[md$Replay_type=="aversive"]
md$Ysig[md$Replay_type=="differential"] <- md$Ysig_rewarding[md$Replay_type=="differential"] - md$Ysig_aversive[md$Replay_type=="differential"]

md$Mmult <- rep(NA,nrow(md))
md$Mmult[md$Replay_type=="rewarding"] <- md$Mmult_rewarding[md$Replay_type=="rewarding"]
md$Mmult[md$Replay_type=="aversive"] <- md$Mmult_aversive[md$Replay_type=="aversive"]
md$Mmult[md$Replay_type=="differential"] <- md$Mmult_rewarding[md$Replay_type=="differential"] - md$Mmult_aversive[md$Replay_type=="differential"]

md$Msig <- rep(NA,nrow(md))
md$Msig[md$Replay_type=="rewarding"] <- md$Msig_rewarding[md$Replay_type=="rewarding"]
md$Msig[md$Replay_type=="aversive"] <- md$Msig_aversive[md$Replay_type=="aversive"]
md$Msig[md$Replay_type=="differential"] <- md$Msig_rewarding[md$Replay_type=="differential"] - md$Msig_aversive[md$Replay_type=="differential"]


md$PathNum <- rep(NA,nrow(md))
md$PathNum[md$Replay_type=="rewarding" & md$Path1Type=="rewarding"] <- 1
md$PathNum[md$Replay_type=="rewarding" & md$Path1Type=="aversive"] <- 2
md$PathNum[md$Replay_type=="aversive" & md$Path1Type=="rewarding"] <- 2
md$PathNum[md$Replay_type=="aversive" & md$Path1Type=="aversive"] <- 1
md$PathNum <- as.factor(md$PathNum)

md$PathVal <- rep(NA,nrow(md))
md$PathVal[md$Replay_type=="rewarding"] <- apply(select(filter(md,Replay_type=="rewarding"),nV_1:nV_2),1,max)
md$PathVal[md$Replay_type=="aversive"] <- apply(select(filter(md,Replay_type=="aversive"),nV_1:nV_2),1,min)

md$PathProb <- rep(NA,nrow(md))
md$PathProb[md$Replay_type=="rewarding"] <- md$RewProb[md$Replay_type=="rewarding"]
md$PathProb[md$Replay_type=="aversive"] <- round(1-md$RewProb[md$Replay_type=="aversive"],1)

md$PathCertainty <- md$PathProb
md$PathCertainty[which(md$PathProb==0.1)] <- rep(0.9,sum(md$PathProb==0.1,na.rm=TRUE))
md$PathCertainty[which(md$PathProb==0.3)] <- rep(0.7,sum(md$PathProb==0.3,na.rm=TRUE))

md$PathRecency <- rep(NA,nrow(md))
md$PathRecency[md$Replay_type=="rewarding"] <- md$RewPathRecency_block[md$Replay_type=="rewarding"]
md$PathRecency[md$Replay_type=="aversive"] <- md$LossPathRecency_block[md$Replay_type=="aversive"]
md$PathRecency <- log(md$PathRecency+1)

md$Choice <- as.factor(md$Choice)

md <- filter(md,
             Lag>10,Lag<100,
             Replay_type=="rewarding" | Replay_type=="aversive")

md <- md[abs(scale(md$Ysig))<10,]

md$PathVal <- scale(md$PathVal,center=TRUE,scale=FALSE)
md$RT <- scale(md$RT,center=TRUE,scale=FALSE)
md$PathProb <- scale(md$PathProb,center=TRUE,scale=FALSE)
md$PathRecency <- scale(md$PathRecency,center=TRUE,scale=FALSE)
md$Ymult <- scale(md$Ymult,center=TRUE,scale=FALSE)
md$Ysig <- scale(md$Ysig,center=TRUE,scale=FALSE)
md$Sequenceness <- scale(md$Sequenceness,center=TRUE,scale=FALSE)

# MODEL 9
m9 <- lmer(Ysig ~ Choice*Replay_type + RT + (1|Subject),
            data=filter(md,Lag==md$Lag[1]), REML=FALSE,
            control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

summary(m9)
cat_plot(m9,pred=Choice,modx=Replay_type)

e <- emmeans(m9,pairwise~Choice*Replay_type)
e <- as.data.frame(e$emmeans)

# MODEL 10
m10 <- lmer(Sequenceness ~ Replay_type*Choice + Ysig + RT + (1|Subject/Lag),
           data=md, REML=FALSE,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

summary(m10)

jtools::effect_plot(m10,pred=Ysig,interval=TRUE) + 
  theme_classic() + 
  coord_cartesian(xlim=c(-0.02,0.02),ylim=c(-0.03,0.03))

cat_plot(m10,pred=Choice,modx=Replay_type)

e <- emmeans(m10,pairwise~Choice*Replay_type)
e <- as.data.frame(e$emmeans)

# ===========================================================================
# Hybrid models with replay prediction
# ===========================================================================

hd <- read.csv('D:/2020_RiskyReplay/results/modelling/hybridmodelswithreplay.csv')

hd$modelNum <- as.factor(hd$modelNum)
hd$Subject <- as.factor(hd$Subject)

# --------

md <- hd %>%
  filter()

m <- lmer(Replay ~ replayinitiated*BIC*modelNum + (1|Subject/Lag),
          data=hd, REML=FALSE,
          control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))
summary(m)

pd <- hd %>%
  filter(as.integer(modelNum)>4) %>%
  group_by(Subject,Lag,modelNum) %>%
  mutate(corr=cor.test(replayinitiated,Replay,na.rm=TRUE)[4][[1]]) %>%
  group_by(Subject,modelNum) %>%
  summarise(corr=mean(corr),BIC=mean(BIC)) %>%
  arrange(Subject,BIC) %>%
  group_by(Subject) %>%
  mutate(modelOrder=1:n(),
         maxModel=BIC==max(BIC))

fig <- ggplot(pd, aes(x=BIC,y=corr,color=Subject)) + 
  geom_point() + 
  facet_wrap(~modelNum)
ggplotly(fig)

fig <- ggplot(filter(pd,maxModel==TRUE), aes(x=modelNum,y=corr,color=Subject)) + 
  geom_point()
ggplotly(fig)