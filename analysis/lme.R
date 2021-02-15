library(plyr)
library(dplyr)
library(Hmisc)
library(reshape2)
library(e1071)
library(gtools)
library(lmerTest)
library(lsmeans)
library(effects)
library(car)
library(interactions)
library(viridis)
library(ggplot2)
library(svglite)
library(jtools)
library(gridExtra)
# library(ppcor)

findInt <- function(model, value) {
  function(x) {
    predict(model, data.frame(EV=x), type="response") - value
  }
}

# ------------------------------------------------------------------------------

dohalves <- FALSE
avLag <- FALSE
removeOutliers <- TRUE

epoch <- "decision" # "decision", "transition", or "outcome"

# ------------------------------------------------------------------

qs <- read.csv("questionnaire_PCA.csv")

# ------------------------------------------------------------------------------

lagRange <- c(20,80)

setwd("D:/Files/2020_RiskyReplay/2021")

# --- load planning period
if (epoch == "decision" | epoch=="transition"){
  if (!dohalves) {
    d <- read.csv(paste("longSeq_decision_withCatch_withBetas.csv",sep=""))
  } else if (dohalves & epoch=="decision") {
    tmp1 <- read.csv(paste("longSeq_decision-firsthalf.csv",sep=""))
    tmp1$Part <- rep("first",dim(tmp1)[1])
    tmp2 <- read.csv(paste("longSeq_decision-lasthalf.csv",sep=""))
    tmp2$Part <- rep("last",dim(tmp2)[1]) 
    d <- rbind(tmp1,tmp2)
  }
  d <- filter(d,Direction!=2,Lag>lagRange[1],Lag<lagRange[2]) # to reduce file size
}

# --- transition planning period
if (epoch == "transition") {
  trans <- read.csv(paste("longSeq_transition.csv",sep=""))
  trans <- filter(trans,Direction!=2,Lag>lagRange[1],Lag<lagRange[2])
  
  trans_in_d <- match(do.call("paste", select(d,Subject,ExpTrial,Path,Direction,Lag)),
                      do.call("paste", select(trans,Subject,ExpTrial,Path,Direction,Lag)))
  
  d$TransSeq <- rep(NA,dim(d)[1])
  d$Selected <- rep(NA,dim(d)[1])
  d$Expected <- rep(NA,dim(d)[1])
  
  d$TransSeq[!is.na(trans_in_d)] = trans$Seq[trans_in_d[!is.na(trans_in_d)]]
  d$Selected[!is.na(trans_in_d)] = trans$Selected[trans_in_d[!is.na(trans_in_d)]]
  d$Expected[!is.na(trans_in_d)] = trans$Expected[trans_in_d[!is.na(trans_in_d)]]
}

# --- value outcome period
if (epoch == "outcome") {
  d <- read.csv(paste("longSeq_outcome.csv",sep=""))
  d <- filter(d,Direction!=2,Lag>lagRange[1],Lag<lagRange[2]) # to reduce file size
}

# questionnaires
qs$AvRisk <- rowMeans(select(qs,Ethical:Social))

d$RiskSeeking <- rep(NA,dim(d)[1])
d$Anxiety <- rep(NA,dim(d)[1])
d$IUS <- rep(NA,dim(d)[1])
d$Worry <- rep(NA,dim(d)[1])
d$Ethical <- rep(NA,dim(d)[1])
d$Financial <- rep(NA,dim(d)[1])
d$Health <- rep(NA,dim(d)[1])
d$Recreational <- rep(NA,dim(d)[1])
d$Social <- rep(NA,dim(d)[1])
d$AvRisk <- rep(NA,dim(d)[1])
for (s in unique(d$Subject)) {
  d$RiskSeeking[d$Subject==s] = rep(qs$RiskSeeking[qs$Subject==s],sum(d$Subject==s))
  d$Anxiety[d$Subject==s] = rep(qs$Anxiety[qs$Subject==s],sum(d$Subject==s))
  d$IUS[d$Subject==s] = rep(qs$IUS[qs$Subject==s],sum(d$Subject==s))
  d$Worry[d$Subject==s] = rep(qs$Worry[qs$Subject==s],sum(d$Subject==s))
  d$Ethical[d$Subject==s] = rep(qs$Ethical[qs$Subject==s],sum(d$Subject==s))
  d$Financial[d$Subject==s] = rep(qs$Financial[qs$Subject==s],sum(d$Subject==s))
  d$Health[d$Subject==s] = rep(qs$Health[qs$Subject==s],sum(d$Subject==s))
  d$Recreational[d$Subject==s] = rep(qs$Recreational[qs$Subject==s],sum(d$Subject==s))
  d$Social[d$Subject==s] = rep(qs$Social[qs$Subject==s],sum(d$Subject==s))
  d$AvRisk[d$Subject==s] = rep(qs$AvRisk[qs$Subject==s],sum(d$Subject==s))
}

# CLEAN UP
d$Subject = as.factor(as.character(d$Subject))

d$Path[d$Path==1] = "rewarding"
d$Path[d$Path==2] = "aversive"
d$Path[d$Path==3] = "differential"
d$Path = as.factor(d$Path)

d$Choice[d$Choice==1] = "risky"
d$Choice[d$Choice==2] = "safe"
d$Choice = as.factor(d$Choice)

d$Direction[d$Direction==1] = "fwd"
d$Direction[d$Direction==2] = "bwd"
d$Direction[d$Direction==3] = "fwdbwd"
d$Direction = as.factor(d$Direction)

d$Forced <- as.factor(d$Forced != 0)

d$OverallBlock <- d$Block
d$OverallBlock[d$Practice==1] = 0
d$OverallBlock <- d$OverallBlock + 1

d <- d %>%
  group_by(Subject,OverallBlock) %>%
  mutate(BlockAcc=mean(Acc,na.rm=T))

if (avLag) {
  d <- d %>%
    filter(Lag > lagRange[1], Lag < lagRange[2]) %>%
    group_by(Subject,ExpTrial,Path,Direction) %>%
    mutate(Seq = mean(Seq),
           TransSeq = mean(TransSeq)) %>%
    filter(Lag == lagRange[1]+10) %>%
    as.data.frame()
}


# skewness(d$Seq)

#-------------------------------------------
# Add variables
#-------------------------------------------

d$RewProb <- d$P
d$RewProb[d$nV_2>0] = 1-d$P[d$nV_2>0]

d$pEV <- rep(NA,dim(d)[1])
idx <- d$Path=="rewarding" | d$Path=="differential"
d$pEV[idx] = d$RewVal[idx] * d$PathProb[idx]
idx <- d$Path=="aversive"
d$pEV[idx] = d$LossVal[idx] * d$PathProb[idx]

d$Certainty <- d$P
d$Certainty[d$P==.1] = .9
d$Certainty[d$P==.3] = .7

d$RewProb <- round(d$RewProb,1)
d$PathProb <- round(d$PathProb,1)

d$CatRewProb <- round(d$RewProb,1)
d$CatRewProb[d$RewProb==.1 | d$RewProb==.3] = "A_unlikely"
d$CatRewProb[d$RewProb==.5] = "B_equal"
d$CatRewProb[d$RewProb==.7 | d$RewProb==.9] = "C_likely"

d$CatProb <- round(d$PathProb,1)
d$CatProb[d$CatProb==.1 | d$CatProb==.3] = "A_unlikely"
d$CatProb[d$CatProb==.5] = "B_equal"
d$CatProb[d$CatProb==.7 | d$CatProb==.9] = "C_likely"

if (epoch == "transition") {
  d$Selected[which(d$Path=="differential" & d$Selected==1)] = "rewarding"
  d$Selected[which(d$Path=="differential" & d$Selected==2)] = "aversive"
  d$Selected[which(d$Path!="differential" & d$Selected==0)] = FALSE
  d$Selected[which(d$Path!="differential" & d$Selected==1)] = TRUE
}

d$BinChoice <- rep(0,dim(d)[1])
d$BinChoice[d$Choice=="risky"] = 1

d$prev_Choice[d$prev_Choice==1] = "risky"
d$prev_Choice[d$prev_Choice==2] = "safe"
d$prev_Choice <- as.factor(d$prev_Choice)

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

d <- d %>% group_by(Subject) %>% mutate(BlockOrder = getmode(nV_1[Block < 6])>0)
d$BlockOrder[d$BlockOrder==T] = "rew_av"
d$BlockOrder[d$BlockOrder==F] = "av_rew"


d$loseStreak[is.na(d$loseStreak)] = 0
d$winStreak[is.na(d$winStreak)] = 0
d$winloseStreak <- d$winStreak
d$winloseStreak[d$winStreak == 0] = -d$loseStreak[d$winStreak == 0]

d$CatStreak <- rep("avoiding",dim(d)[1])
d$CatStreak[d$winStreak > 0] = "winning"
d$CatStreak[d$loseStreak > 0] = "losing"

d <- d %>% group_by(Subject) %>% mutate(CatRewVal=RewVal>median(RewVal),
                                        CatLossVal=LossVal<median(LossVal)) %>% as.data.frame()
d$CatRewVal[d$CatRewVal==TRUE] = "Large"
d$CatRewVal[d$CatRewVal==FALSE] = "Small"
d$CatLossVal[d$CatLossVal==TRUE] = "Large"
d$CatLossVal[d$CatLossVal==FALSE] = "Small"

d$CatValDiff <- rep(NA,dim(d)[1])
d$CatValDiff[d$CatRewVal=="Small" & d$CatLossVal=="Small"] = "sR_sL"
d$CatValDiff[d$CatRewVal=="Small" & d$CatLossVal=="Large"] = "sR_lL"
d$CatValDiff[d$CatRewVal=="Large" & d$CatLossVal=="Small"] = "lR_sL"
d$CatValDiff[d$CatRewVal=="Large" & d$CatLossVal=="Large"] = "lR_lL"


d$rEV <- d$RewVal*d$RewProb
d$lEV <- d$LossVal*(1-d$RewProb)

d$ActualPath <- rep(NA,dim(d)[1])
d$ActualPath[d$Path=="rewarding" & d$nV_1 > d$nV_2] = 0 # path 1
d$ActualPath[d$Path=="rewarding" & d$nV_2 > d$nV_1] = 1 # path 2
d$ActualPath[d$Path=="aversive"  & d$nV_1 < d$nV_2] = 0
d$ActualPath[d$Path=="aversive"  & d$nV_2 < d$nV_1] = 1

d$Val <- rep(NA,dim(d)[1])
d$Val[d$ActualPath==0 & d$Path!="differential"] = d$nV_1[d$ActualPath==0 & d$Path!="differential"]
d$Val[d$ActualPath==1 & d$Path!="differential"] = d$nV_2[d$ActualPath==1 & d$Path!="differential"]

d$rewSeq <- rep(NA,dim(d)[1])
d$rewSeq[d$Path=="differential"] = d$Seq[d$Path=="rewarding"]
d$lossSeq <- rep(NA,dim(d)[1])
d$lossSeq[d$Path=="differential"] = d$Seq[d$Path=="aversive"]

if (epoch=="transition"){
  d$TransRewSeq <- rep(NA,dim(d)[1])
  d$TransRewSeq[d$Path=="differential"] = d$TransSeq[d$Path=="rewarding"]
  d$TransLossSeq <- rep(NA,dim(d)[1])
  d$TransLossSeq[d$Path=="differential"] = d$TransSeq[d$Path=="aversive"]
}

d$CatChoice <- rep(NA,dim(d)[1])
d$CatChoice[d$Choice=="risky" & d$Acc==1] = "hit"
d$CatChoice[d$Choice=="risky" & d$Acc==0] = "FA"
d$CatChoice[d$Choice=="safe" & d$Acc==1] = "CR"
d$CatChoice[d$Choice=="safe" & d$Acc==0] = "miss"
d$CatChoice <- as.factor(d$CatChoice)

d$BlockEV_risky <- rep(NA,dim(d)[1])
d$BlockEV_safe <- rep(NA,dim(d)[1])
for (s in unique(d$Subject)) {
  for (block in unique(d$OverallBlock)) {
    d$BlockEV_risky[d$Subject==s & d$OverallBlock==block] = mean(d$EV[d$Subject==s & d$OverallBlock==block & d$Choice=="risky"])
    d$BlockEV_safe[d$Subject==s & d$OverallBlock==block] = mean(d$EV[d$Subject==s & d$OverallBlock==block & d$Choice=="safe"])
  }
}

if (epoch == "transition"){
  d$TransitionType <- rep(NA,dim(d)[1])
  d$TransitionType[(d$Transition==1 & d$nV_1 > 0) | 
                     (d$Transition==2 & d$nV_2 > 0)] = "rewarding"
  d$TransitionType[(d$Transition==1 & d$nV_1 < 0) | 
                     (d$Transition==2 & d$nV_2 < 0)] = "aversive"
  d$TransitionType <- as.factor(d$TransitionType)
}
  
d$negChange <- rep(NA,dim(d)[1])
d$probChange <- rep(NA,dim(d)[1])
for (s in unique(d$Subject)) {
  for (block in unique(d$OverallBlock)) {
    thisd <- filter(d,
                    Subject==s,
                    OverallBlock==block,
                    Lag==unique(d$Lag)[1],
                    Path==unique(d$Path[1]),
                    Direction==unique(d$Direction)[1])
    trls <- thisd$ExpTrial
    for (trl in 2:length(trls)) {
      if (trls[trl] - trls[trl-1] == 1) {
        negChange <- thisd$nCombo[trl] != thisd$nCombo[trl-1]
        probChange <- thisd$P[trl] != thisd$P[trl-1]
        idx <- d$Subject == s &
          d$OverallBlock == block &
          d$ExpTrial == trls[trl]
        d$negChange[idx] = negChange
        d$probChange[idx] = probChange
      }
    }
  }
}

# d <- d %>%
#   group_by(Subject,Direction,Path,Lag) %>%
#   mutate(Outlier = abs(scale(AvSeqAvLag)) > 3)

if (epoch == "transition") {
  d$AvTransSeq <- rep(NA,dim(d)[1])
  for (s in unique(d$Subject)) {
    expTrial <- unique(d$ExpTrial[d$Subject==s])
    print(s)
    for (t in expTrial) {
      for (di in unique(d$Direction)) {
        for (l in unique(d$Lag)) {
          dval <- mean(d$TransSeq[d$Path!="differential" & d$Subject==s & d$ExpTrial==t & d$Direction==di & d$Lag==l])
          d$AvTransSeq[d$Subject==s & d$ExpTrial==t & d$Direction==di & d$Lag==l] = dval
        }
      }
    }
  }
}

# d <- d %>%
#   group_by(Subject,Direction,Path,Lag) %>%
#   mutate(TransOutlier = abs(scale(AvTransSeq)) > 3)


d$AccDiff <- d$EV
d$AccDiff[d$Choice == "risky" & d$Acc==1 & !is.na(d$Acc)]  = d$EV[d$Choice == "risky" & d$Acc==1 & !is.na(d$Acc)] - 1
d$AccDiff[d$Choice == "risky" & d$Acc==0 & !is.na(d$Acc)]  = d$EV[d$Choice == "risky" & d$Acc==0 & !is.na(d$Acc)] + 1
d$AccDiff[d$Choice == "safe" & d$Acc==1 & !is.na(d$Acc)] = (d$EV[d$Choice == "safe" & d$Acc==1 & !is.na(d$Acc)] * (-1)) + 1
d$AccDiff[d$Choice == "safe" & d$Acc==0 & !is.na(d$Acc)] = (d$EV[d$Choice == "safe" & d$Acc==0 & !is.na(d$Acc)] - 1) * (-1)

d$IP <- rep(NA,dim(d)[1])
d$bIP <- rep(NA,dim(d)[1])
subjects <- unique(d$Subject)
for (s in subjects) {
  
  thisd <- filter(d,
                  Subject==s,
                  Forced==F,
                  Path==unique(d$Path)[1],
                  Direction==unique(d$Direction)[1],
                  Lag==unique(d$Lag)[1])
  
  L <- glm(BinChoice ~ EV, family="binomial", data=thisd)
  thisIP <- uniroot(findInt(L, .5), range(seq(-100,100,0.1)))$root
  d$IP[d$Subject==s] = rep(thisIP,sum(d$Subject==s))
  
  for (b in unique(thisd$OverallBlock)) {
    bd <- filter(thisd,OverallBlock==b)
    L <- glm(BinChoice ~ EV, family="binomial", data=bd)
    thisIP <- tryCatch({
      uniroot(findInt(L, .5), range(seq(min(d$EV),max(d$EV),0.1)))$root
    }, error=function(cond){
      return(NA)
    })
    d$bIP[d$Subject==s & d$OverallBlock==b] = rep(thisIP,sum(d$Subject==s & d$OverallBlock==b))
  }
}

d$pType <- rep(NA,dim(d)[1])
d$pType[d$nV_1 > 0 & d$nV_2 < 0] = 1
d$pType[d$nV_1 < 0 & d$nV_2 > 0] = 2
d$pType[d$nV_1 > 0 & d$nV_2 > 0] = 3
d$pType[d$nV_1 < 0 & d$nV_2 < 0] = 4
d$pType <- as.factor(as.character(d$pType))

d$DomBlock <- rep(NA,dim(d)[1])
d$DomBlock[d$OverallBlock < 7 & d$BlockOrder == "rew_av"] = "rewarding"
d$DomBlock[d$OverallBlock > 6 & d$BlockOrder == "rew_av"] = "aversive"
d$DomBlock[d$OverallBlock < 7 & d$BlockOrder == "av_rew"] = "aversive"
d$DomBlock[d$OverallBlock > 6 & d$BlockOrder == "av_rew"] = "rewarding"


# Add experience count (overall & per block)

d$prev_Outcome[is.na(d$prev_Outcome)] = 0

d <- d %>%
  group_by(Subject,Direction,Path,Lag) %>%
  mutate(p1_exp = cumsum(Transition==1),
         p2_exp = cumsum(Transition==2),
         rew_exp = cumsum(prev_Outcome > 0 & Choice=="risky"),
         loss_exp = cumsum(prev_Outcome < 0 & Choice=="risky"))

d <- d %>%
  group_by(Subject,Direction,Path,Lag,OverallBlock) %>%
  mutate(p1_exp_block = cumsum(Transition==1),
         p2_exp_block = cumsum(Transition==2),
         rew_exp_block = cumsum(prev_Outcome > 0 & Choice=="risky"),
         loss_exp_block = cumsum(prev_Outcome < 0 & Choice=="risky"))

tmp1 <- d %>%
  group_by(Subject,Direction,Path,Lag) %>%
  filter(OverallBlock < 7) %>%
  mutate(p1_exp_half = cumsum(Transition==1),
         p2_exp_half = cumsum(Transition==2),
         rew_exp_half = cumsum(prev_Outcome > 0 & Choice=="risky"),
         loss_exp_half = cumsum(prev_Outcome < 0 & Choice=="risky"))

tmp2 <- d %>%
  group_by(Subject,Direction,Path,Lag) %>%
  filter(OverallBlock > 6) %>%
  mutate(p1_exp_half = cumsum(Transition==1),
         p2_exp_half = cumsum(Transition==2),
         rew_exp_half = cumsum(prev_Outcome > 0 & Choice=="risky"),
         loss_exp_half = cumsum(prev_Outcome < 0 & Choice=="risky"))

d <- rbind(tmp1,tmp2)

#-------------------------------------------
# Plot behaviour
#-------------------------------------------

pd <- d %>% 
  filter(Forced==FALSE,
         Subject!="263098",
         Subject!="680913",
         Lag==unique(d$Lag)[1],
         Path==unique(d$Path[1]))

rd <- d %>% 
  filter(Forced==FALSE,
         Subject=="263098" | Subject=="680913",
         Lag==unique(d$Lag)[1],
         Path==unique(d$Path[1]))

for (i in c(1,2,3)) {
  
  if (i == 1){
    p <- pd %>% group_by(Subject,Choice) %>% summarise(mean=mean(Acc))
    r <- rd %>% group_by(Subject,Choice) %>% summarise(mean=mean(Acc))
  } else if (i == 2){
    p <- pd %>% group_by(Subject,Choice) %>% summarise(mean=mean(EV))
    r <- rd %>% group_by(Subject,Choice) %>% summarise(mean=mean(EV))
  } else if (i == 3){
    p <- pd %>% group_by(Subject,Choice) %>% summarise(mean=mean(RT))
    r <- rd %>% group_by(Subject,Choice) %>% summarise(mean=mean(RT))
  }
  dc <- dcast(p,Subject~Choice,mean)
  t.test(dc$risky,dc$safe,paired=T)
  
  ggplot() + 
    geom_point(data=p, aes(x=Choice,y=mean,group=Subject,color=Subject),size=2) + 
    geom_line(data=p, aes(x=Choice,y=mean,group=Subject,color=Subject)) + 
    geom_point(data=r, aes(x=Choice,y=mean,group=Subject),size=2,colour="grey") + 
    geom_line(data=r, aes(x=Choice,y=mean,group=Subject),alpha=0.5,colour="grey") + 
    theme_classic() + scale_color_viridis(discrete=T)
  
}

#-------------------------------------------
# BEHAVIOURAL MODELS
#-------------------------------------------

# acc & rt graphs
md <- filter(d,Direction=="fwdbwd",
             Lag == unique(d$Lag)[1],
             Forced==FALSE,
             Path==unique(d$Path)[1]
)

dtable <- md %>% 
  group_by(Subject,Choice) %>% 
  summarise(
    acc = mean(Acc),
    ev = mean(EV),
    rt = mean(RT)
  ) %>% 
  as.data.frame()

ggplot(data=dtable, aes(x=Choice,y=acc,group=Subject,color=Subject,fill=Subject)) + 
  geom_point(size=3) + 
  geom_line() +
  theme_classic() + 
  scale_color_viridis(discrete=T) + scale_fill_viridis(discrete=T) + ylim(.4,1)

ggplot(data=dtable, aes(x=Choice,y=ev,group=Subject,color=Subject,fill=Subject)) + 
  geom_point(size=3) + 
  geom_line() +
  theme_classic() + 
  scale_color_viridis(discrete=T) + scale_fill_viridis(discrete=T) + ylim(-3,3)

ggplot(data=dtable, aes(x=Choice,y=rt,group=Subject,color=Subject,fill=Subject)) + 
  geom_point(size=3) + 
  geom_line() +
  theme_classic() + 
  scale_color_viridis(discrete=T) + scale_fill_viridis(discrete=T) + ylim(5,16)

# modelling

md <- filter(d,Direction=="fwdbwd",
             Lag == unique(d$Lag)[1],
             Forced==FALSE,
             Path=="differential",
             Subject!="263098",
             Subject!="680913"
)
md$RewVal <- scale(md$RewVal,scale=F,center=T)
md$LossVal <- scale(abs(md$LossVal),scale=F,center=T)
md$RewProb <- scale(md$RewProb,scale=F,center=T)
md$RiskSeeking <- scale(md$RiskSeeking,scale=F,center=T)
md$Anxiety <- scale(md$Anxiety,scale=F,center=T)

m1 <- glmer(BinChoice ~ (1|Subject) + 
              RewVal*LossVal*RewProb, 
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

summary(m1)
vif(m1)


interact_plot(m1,pred=RewProb,modx=RewVal,
              interval=T,modx.values=quantile(md$RewVal,c(.05,.5,.95)),
              vary.lty=F,colors=c("#C970FF","#00B2FF","#1EFFA0")) + 
  theme_classic() + ylim(0,1)
interact_plot(m1,pred=RewProb,modx=LossVal,
              interval=T,modx.values=quantile(md$LossVal,c(.05,.5,.95)),
              vary.lty=F,colors=c("#FFD500","#FF7C00","#FF0087")) + 
  theme_classic() + ylim(0,1)

ss <- sim_slopes(m1,pred=rEV,modx=lEV,data=md,
                 modx.values=quantile(md$lEV,c(.25,.5,.75)))


m2 <- glmer(BinChoice ~ (RT|Subject) + 
              RewVal*LossVal*RewProb*RT, 
            data=md, family="binomial",
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

summary(m2)

# m2 <- glmer(BinChoice ~ (1|Subject) + 
#               EV*Anxiety*Certainty + EV*RiskSeeking*Certainty, 
#             data=md, family="binomial",
#             control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))
# 
# summary(m2)
# vif(m2)
# 
# interact_plot(m2,pred=EV,modx=Anxiety,modx.values=unique(md$Anxiety),
#               interval=F,vary.lty=F) + theme_classic() + scale_color_viridis(option="C")

# indifference points

S <- unique(md$Subject)
N <- length(S)
evrange <- seq(-100,100,0.1)
IP <- data.frame(Subject=S,
                 Point = rep(NA,N),
                 pval = rep(NA,N))

for (subject in S) {
  
  thisd <- md[md$Subject==subject,]
  L <- glm(Choice ~ EV, family="binomial", data=thisd)
  
  thisIP = uniroot(findInt(L, .5), range(evrange))$root
  IP$Point[IP$Subject==subject] <- thisIP
  
  IP$pval[IP$Subject==subject] <- coef(summary(L))[2,4]
  
  pred <- predict(L, data.frame(EV=evrange), type="response", se.fit=TRUE)
  P <- data.frame(Subject = rep(subject,length(evrange)),
                  EV = evrange,
                  pred = pred$fit,
                  lower = pred$fit - (2*pred$se.fit),
                  upper = pred$fit + (2*pred$se.fit))
  if (subject==S[1]){
    plotfit <- P
  } else {
    plotfit <- rbind(plotfit,P)
  }
}

plotfit <- plotfit[order(plotfit$Subject,plotfit$EV),]

thisd <- md[md$Subject!="263098" & d$Subject!="680913",]
L <- glm(Choice ~ EV, family="binomial", data=thisd)
grandPred <- predict(L, data.frame(EV=evrange), type="response", se.fit=TRUE)

plotfit <- rbind(plotfit,
                 data.frame(Subject=rep("grand",length(evrange)),
                            EV = evrange,
                            pred = grandPred$fit,
                            lower = grandPred$fit - (2*grandPred$se.fit),
                            upper = grandPred$fit + (2*grandPred$se.fit)))

# IP <- IP[IP$pval < .05,]

ggplot() + 
  geom_line(data=plotfit,
            aes(x=EV,y=lower,group=Subject),
            linetype="dashed") +
  geom_line(data=plotfit,
            aes(x=EV,y=upper,group=Subject),
            linetype="dashed") +
  geom_line(data=plotfit,
            aes(x=EV,y=pred,group=Subject)) + 
  geom_point(data=IP,
             aes(x=Point,y=0.5,group=Subject),
             size=3) + 
  geom_hline(yintercept=0.5, alpha=0.2) + 
  geom_vline(xintercept=1, alpha=0.2) + 
  xlim(-9,11) + ylim(0,1) + 
  facet_wrap(~Subject,nrow=3,ncol=5) + 
  theme_classic()


#-------------------------------------------
# REPLAY-ONLY MODELS
#-------------------------------------------

# (VAL & PROB)
md <- filter(d,Direction=="fwdbwd",
             Lag > lagRange[1], Lag < lagRange[2],
             Forced==F,
             Path!="differential",
             Subject!="263098",
             Subject!="680913"
)

md$Path <- droplevels(md$Path)
md$Subject <- droplevels(md$Subject)

md <- md %>%
  group_by(Subject) %>%
  mutate(Seq = scale(Seq)) %>% # standardize sequenceness per participant
  group_by(Subject,ExpTrial) %>%
  mutate(Outlier = max(abs(Seq)) > 3) # remove outliers using average sequenceness

outliers <- md %>% 
  group_by(Subject) %>% 
  filter(Outlier==T,Lag==unique(md$Lag)[1],Path==unique(md$Path)[1]) %>% 
  summarise(n = n(), small=sum(Seq < -3), big=sum(Seq > 3))

if (removeOutliers) {
  md <- filter(md,Outlier==F)
}

md$Val <- scale(md$Val,center=T,scale=F)
md$sRewProb <- scale(md$RewProb,center=T,scale=F)
md$PathProb <- scale(md$PathProb,center=T,scale=F)
md$ExpTrial <- scale(md$ExpTrial,center=T,scale=F)
md$Lag <- scale(md$Lag,center=T,scale=F)
md$EV <- scale(md$EV,center=T,scale=F)
md$pEV <- scale(md$pEV,center=T,scale=F)
md$Seq <- scale(md$Seq,center=T,scale=F)


m1 <- lmer(Seq ~ (1|Subject) + (1|Subject:ExpTrial) + (1|Subject:ExpTrial:Lag) +
             Val*PathProb, 
           data=md, REML=F,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

summary(m1)

round(select(as.data.frame(summary(m1)$coefficients),Estimate,`Pr(>|t|)`),3)

interact_plot(m1,pred=Val,
              modx=PathProb,modx.values=quantile(md$PathProb,seq(0,1,.2)),
              interval=T,vary.lty=F) + theme_classic() + scale_color_viridis() + scale_fill_viridis() + 
  ylim(-0.25,0.25)



m2 <- lmer(Seq ~ (1|Subject) + (1|Subject:ExpTrial) + (1|Subject:ExpTrial:Lag) +
             Val*Choice, 
           data=md, REML=F,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

summary(m2)

interact_plot(m2,pred=Val,modx=Choice,interval=T) + theme_classic() + ylim(-0.2,0.2)

ss <- sim_slopes(m2,pred=Val,modx=Choice)




#----------------------------------------
# (CHOICE AS DV)

md <- filter(d,Direction=="fwdbwd",
             Lag > lagRange[1], Lag < lagRange[2],
             Forced==FALSE,
             Path=="differential",
             Subject!="263098",
             Subject!="680913"
)

md <- md %>%
  group_by(Subject) %>%
  mutate(Seq = scale(Seq),
         AvSeq = scale(AvSeq)) %>% # standardize sequenceness per participant
  group_by(Subject,ExpTrial) %>%
  mutate(Outlier = max(abs(Seq)) > 3 | max(abs(AvSeq)) > 3) # remove outliers using average sequenceness

outliers <- md %>% 
  group_by(Subject) %>% 
  filter(Outlier==T,Lag==unique(md$Lag)[1],Path==unique(md$Path)[1]) %>% 
  summarise(n = n(), small=sum(AvSeq < -3 | Seq < -3), big=sum(AvSeq > 3 | Seq > 3))

if (removeOutliers) {
  md <- filter(md,Outlier==F)
}

md$Path <- droplevels(md$Path)
md$Subject <- droplevels(md$Subject)

md$LossVal <- scale(abs(md$LossVal),center=T,scale=F)
md$RewVal <- scale(md$RewVal,center=T,scale=F)
md$RewProb <- scale(md$RewProb,center=T,scale=F)
md$PathProb <- scale(md$PathProb,center=T,scale=F)
md$EV <- scale(md$EV,center=T,scale=F)

md$ExpTrial <- scale(md$ExpTrial,center=T,scale=F)
md$OverallBlock <- scale(md$OverallBlock,center=T,scale=F)

md$rewloss_exp <- scale(md$rew_exp - md$loss_exp,center=T,scale=F)
md$rew_exp <- scale(md$rew_exp,center=T,scale=F)
md$loss_exp <- scale(md$loss_exp,center=T,scale=F)

md$rewloss_exp_block <- scale(md$rew_exp_block - md$loss_exp_block,center=T,scale=F)
md$rew_exp_block <- scale(md$rew_exp_block,center=T,scale=F)
md$loss_exp_block <- scale(md$loss_exp_block,center=T,scale=F)

md$rewloss_exp_half <- scale(md$rew_exp_half - md$loss_exp_half,center=T,scale=F)
md$rew_exp_half <- scale(md$rew_exp_half,center=T,scale=F)
md$loss_exp_half <- scale(md$loss_exp_half,center=T,scale=F)

md$Seq <- scale(md$Seq,center=T,scale=F)
md$AvSeq <- scale(md$AvSeq,center=T,scale=F)


m1 <- glmer(BinChoice ~ (1|Subject/Lag) +
              Seq + AvSeq, 
            data=md, family=binomial,
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

summary(m1)
round(select(as.data.frame(summary(m1)$coefficients),Estimate,`Pr(>|z|)`),3)

m2 <- glmer(BinChoice ~ (1|Subject/Lag) +
              Seq*RewVal*LossVal*RewProb + AvSeq*RewVal*LossVal*RewProb, 
            data=md, family=binomial,
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

summary(m2)
round(select(as.data.frame(summary(m2)$coefficients),Estimate,`Pr(>|z|)`),3)

effect_plot(m2,pred=Seq,interval=T) + theme_classic() + ylim(.4,.8) + xlim(-3,3)


interact_plot(m2,pred=RewProb,
              modx=RewVal,modx.values=quantile(md$RewVal,c(.05,.5,.95)),
              mod2=Seq,mod2.values=quantile(md$Seq,c(.05,.5,.95)),
              interval=T,vary.lty=F,colors=c("#C970FF","#00B2FF","#1EFFA0")) + 
  theme_classic()

interact_plot(m2,pred=RewProb,
              modx=LossVal,modx.values=quantile(md$LossVal,c(.05,.5,.95)),
              mod2=Seq,mod2.values=quantile(md$Seq,c(.05,.5,.95)),
              interval=T,vary.lty=F,colors=c("#FFD500","#FF7C00","#FF0087")) + 
  theme_classic()


m3 <- glmer(BinChoice ~ (1|Subject/Lag) +
              Seq*EV + AvSeq*EV, 
            data=md, family=binomial,
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))

summary(m3)
round(select(as.data.frame(summary(m3)$coefficients),Estimate,`Pr(>|z|)`),3)

interact_plot(m3,pred=EV,
              modx=AvSeq,modx.values=quantile(md$AvSeq,c(.05,.95)),
              interval=T,vary.lty=F) + theme_classic()

# --- contrasts
# https://stats.stackexchange.com/questions/483488/comparing-simple-slopes-with-one-another-in-r

cr <- pairs(emtrends(m2,
                     ~ Seq | RewVal, var="RewProb",
                     at = data.frame(
                       Seq=quantile(md$Seq,c(.05,.5,.95)),
                       RewVal=quantile(md$RewVal,c(.05,.5,.95))
                     )),simple="each")$`simple contrasts for Seq`

cl <- pairs(emtrends(m2,
                     ~ Seq | LossVal, var="RewProb",
                     at = data.frame(
                       Seq=quantile(md$Seq,c(.05,.5,.95)),
                       LossVal=quantile(md$LossVal,c(.05,.5,.95))
                     )),simple="each")$`simple contrasts for Seq`


m4 <- glmer(BinChoice ~ (1|Subject/Lag) +
              Seq*rew_exp_block + Seq*loss_exp_block + AvSeq*rew_exp_block + AvSeq*loss_exp_block, 
            data=md, family=binomial,
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))
summary(m4)

interact_plot(m4,pred=loss_exp_block,modx=Seq,modx.values=quantile(md$Seq,c(.05,.5,.95)),
              interval=T,vary.lty=F,colors=c("#FF0074","#00F08A")) + theme_classic() + ylim(0.3,1)

interact_plot(m4,pred=rew_exp_block,modx=Seq,modx.values=quantile(md$Seq,c(.05,.5,.95)),
              interval=T,vary.lty=F,colors=c("#FF0074","#00F08A")) + theme_classic() + ylim(0.3,1)

# individual differences

md <- filter(d,Direction=="fwdbwd",
             Lag > lagRange[1], Lag < lagRange[2],
             Forced==FALSE,
             Path=="differential",
             Subject!="263098",
             Subject!="680913"
)

md <- md %>%
  group_by(Subject) %>%
  mutate(zSeq = scale(Seq),
         zAvSeq = scale(AvSeq)) %>% # standardize sequenceness per participant
  group_by(Subject,ExpTrial) %>%
  mutate(Outlier = max(abs(zSeq)) > 3 | max(abs(zAvSeq)) > 3) # remove outliers using average sequenceness

outliers <- md %>% 
  group_by(Subject) %>% 
  filter(Outlier==T,Lag==unique(md$Lag)[1],Path==unique(md$Path)[1]) %>% 
  summarise(n = n(), small=sum(AvSeq < 3), big=sum(AvSeq>3))

if (removeOutliers) {
  md <- filter(md,Outlier==F)
}

md$Path <- droplevels(md$Path)
md$Subject <- droplevels(md$Subject)
md <- md %>%
  group_by(Subject) %>%
  mutate(medRT = median(RT),
         avRT = mean(RT),
         avAcc = mean(Acc),
         avbAcc = mean(bAcc)) %>%
  as.data.frame()

md$RewProb <- scale(md$RewProb,center=T,scale=F)
md$PathProb <- scale(md$PathProb,center=T,scale=F)
md$RiskSeeking = scale(md$RiskSeeking,center=T,scale=F)
md$AvRisk = scale(md$RiskSeeking,center=T,scale=F)
md$Anxiety = scale(md$Anxiety,center=T,scale=F)
md$IP = scale(md$IP,center=T,scale=F)
md$avAcc = scale(md$avAcc,center=T,scale=F)
md$avbAcc = scale(md$avbAcc,center=T,scale=F)
md$avRT = scale(md$avRT,center=T,scale=F)

md$origSeq <- md$Seq
md$Seq <- scale(md$Seq,center=T,scale=F)
md$AvSeq <- scale(md$AvSeq,center=T,scale=F)


# as fixed effects

m0 <- glmer(BinChoice ~ (1|Subject/Lag) +
              EV, 
            data=md, family=binomial,
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))
effect_plot(m0,pred=EV,interval=T) + theme_classic() + xlim(-10,10)

m1 <- glmer(BinChoice ~ (1|Subject/Lag) +
              Seq*RiskSeeking*EV + Seq*Anxiety*EV +
              AvSeq*RiskSeeking*EV + AvSeq*Anxiety*EV, 
            data=md, family=binomial,
            control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))
summary(m1)
round(select(as.data.frame(summary(m1)$coefficients),Estimate,`Pr(>|z|)`),3)
 

# (risk seeking X EV)
interact_plot(m1,pred=EV,
              modx=RiskSeeking,modx.values=unique(md$RiskSeeking),
              interval=F,vary.lty=F) + theme_classic() + scale_color_viridis(option="C")

# (risk seeking X Seq X EV)
interact_plot(m1,pred=EV,
              mod2=RiskSeeking,mod2.values=quantile(md$RiskSeeking,c(.05,.95)),
              modx=Seq,modx.values=quantile(md$Seq,c(.05,.5,.95)),
              interval=T,vary.lty=F,colors=c("#FF005D","#0CF096")) + theme_classic() + ylim(0,1)
# (anxiety X EV)
interact_plot(m1,pred=EV,
              modx=Anxiety,modx.values=unique(md$Anxiety),
              interval=F,vary.lty=F) + theme_classic() + scale_color_viridis(option="C")

# (anxiety X Seq X EV)
interact_plot(m1,pred=EV,
              mod2=Anxiety,mod2.values=quantile(md$Anxiety,c(.05,.95)),
              modx=Seq,modx.values=quantile(md$Seq,c(.05,.5,.95)),
              interval=T,vary.lty=F,colors=c("#FF005D","#0CF096")) + theme_classic() + ylim(0,1)


c_rs <- pairs(emtrends(m1,
                     ~ Seq | RiskSeeking, var="EV",
                     at = data.frame(
                       Seq=quantile(md$Seq,c(.05,.95)),
                       RiskSeeking=quantile(md$RiskSeeking,c(.05,.95))
                     )),simple="each")

c_anx <- pairs(emtrends(m1,
                       ~ Seq | Anxiety, var="EV",
                       at = data.frame(
                         Seq=quantile(md$Seq,c(.05,.95)),
                         Anxiety=quantile(md$Anxiety,c(.25,.75))
                       )),simple="each")


pd <- md %>% group_by(Subject,Choice) %>% summarise(RiskSeeking = mean(RiskSeeking),
                                                       Seq = mean(Seq))
pd <- pd %>% group_by(Subject) %>% mutate(Seq_choiceDiff = lag(Seq)-Seq)
pd <- pd[!is.na(pd$Seq_choiceDiff),]

cor.test(pd$RiskSeeking,pd$Seq_choiceDiff,method="spearman")

ggplot(pd,aes(x=RiskSeeking,y=Seq_choiceDiff,color=RiskSeeking,fill=RiskSeeking)) + 
  geom_point() + 
  geom_smooth(method="lm",se=F,color="black") + 
  theme_classic() + scale_color_viridis(option="C") + scale_fill_viridis(option="C")


pd <- md %>% group_by(Subject,Choice) %>% summarise(RiskSeeking = mean(RiskSeeking),
                                                    Anxiety = mean(Anxiety),
                                                    Seq = mean(Seq))

ggplot(pd,aes(x=Anxiety,y=Seq,color=Anxiety,fill=Anxiety)) + 
  geom_point() + 
  geom_smooth(method="lm",se=F,color="black") + 
  theme_classic() + facet_wrap(~Choice) +
  scale_color_viridis(option="C") + scale_fill_viridis(option="C")


#############################################################################################
# #### TRANSITION

md <- filter(d,Direction=="fwdbwd",
             Lag > lagRange[1], Lag < lagRange[2],
             Choice=="risky",
             Forced==F,
             Path!="differential",
             Subject!="263098",
             Subject!="680913"
)
md <- md[!is.na(md$TransSeq),]

md <- md %>%
  group_by(Subject) %>%
  mutate(TransSeq = scale(TransSeq),
         AvTransSeq = scale(AvTransSeq)) %>% # standardize sequenceness per participant
  group_by(Subject,ExpTrial) %>%
  mutate(Outlier = max(abs(TransSeq)) > 3 | max(abs(AvTransSeq)) > 3) # remove outliers using average sequenceness

outliers <- md %>% 
  group_by(Subject) %>% 
  filter(Outlier==T,Lag==unique(md$Lag)[1],Path==unique(md$Path)[1]) %>% 
  summarise(n = n(), small=sum(AvTransSeq < -3 | TransSeq < -3), big=sum(AvTransSeq > 3 | TransSeq > 3))

if (removeOutliers) {
  md <- filter(md,Outlier==F)
}


md$Path <- droplevels(md$Path)
md$Subject <- droplevels(md$Subject)

md$CatExpec <- md$RewProb > .4

md$RewVal <- scale(md$RewVal,center=T,scale=F)
md$LossVal <- scale(abs(md$LossVal),center=T,scale=F)
md$Expected <- scale(md$Expected,center=T,scale=F)
md$RewProb <- scale(md$RewProb,center=T,scale=F)
md$TransSeq <- scale(md$TransSeq,center=T,scale=F)
md$AvTransSeq <- scale(md$AvTransSeq,center=T,scale=F)
md$RiskSeeking <- scale(md$RiskSeeking,center=T,scale=F)
md$Anxiety <- scale(md$Anxiety,center=T,scale=F)


m1 <- lmer(TransSeq ~ (1|Subject) + (1|Subject:ExpTrial) + (1|Subject:ExpTrial:Lag) +
             Path*TransitionType,
           data=md, REML=F,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))
summary(m1)

cat_plot(m1,pred=TransitionType,modx=Path,
         interval=T,colors=c("#FF005D","#0CF096")) + 
  theme_classic() + ylim(-0.15,.1)

lsmeans(m1,pairwise~Path*TransitionType)


m2 <- lmer(TransSeq ~ (1|Subject) + (1|Subject:ExpTrial) + (1|Subject:ExpTrial:Lag) +
             Path*TransitionType*RiskSeeking + Path*TransitionType*Anxiety,
           data=md, REML=F,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))
summary(m2)

cat_plot(m2,pred=TransitionType,modx=Path,
         interval=T,colors=c("#FF005D","#0CF096")) + 
  theme_classic() + ylim(-0.15,.1)

lsmeans(m2,pairwise~Path*TransitionType)



cat_plot(m2,pred=TransitionType,modx=Path,
         mod2=RiskSeeking,mod2.values=quantile(md$RiskSeeking,c(.05,.95)),
         interval=T,colors=c("#FF005D","#0CF096")) + 
  theme_classic() + ylim(-.4,.3)

cat_plot(m2,pred=TransitionType,modx=Path,
         mod2=Anxiety,mod2.values=quantile(md$Anxiety,c(.05,.95)),
         interval=T,colors=c("#FF005D","#0CF096")) + 
  theme_classic() + ylim(-.4,.3)


interact_plot(m2,pred=Anxiety,modx=Path,
         mod2=TransitionType,
         interval=T,colors=c("#FF005D","#0CF096")) + 
  theme_classic() + ylim(-.4,.3)

c_anx <- pairs(emtrends(m2,
                       ~ TransitionType | Path, var="Anxiety",
                       at = data.frame(
                         Path=c("aversive","rewarding"),
                         TransitionType=c("aversive","rewarding"),
                         Anxiety=quantile(md$Anxiety,c(.05,.95))
                       )),simple="each")


#############################################################################################
# #### OUTCOME

md <- filter(d,Direction=="fwdbwd",
             Lag > lagRange[1], Lag < lagRange[2],
             Choice=="safe",
             Forced==F,
             Path!="differential",
             Subject!="263098",
             Subject!="680913"
)

md <- md %>%
  group_by(Subject) %>%
  mutate(Seq = scale(Seq),
         AvSeq = scale(AvSeq)) %>% # standardize sequenceness per participant
  group_by(Subject,ExpTrial) %>%
  mutate(Outlier = max(abs(Seq)) > 3 | max(abs(AvSeq)) > 3) # remove outliers using average sequenceness

outliers <- md %>% 
  group_by(Subject) %>% 
  filter(Outlier==T,Lag==unique(md$Lag)[1],Path==unique(md$Path)[1]) %>% 
  summarise(n = n(), small=sum(AvSeq < 3), big=sum(AvSeq>3))

if (removeOutliers) {
  md <- filter(md,Outlier==F)
}

m1 <- lmer(Seq ~ (1|Subject) + (1|Subject:ExpTrial) + (1|Subject:ExpTrial:Lag) +
             Path,
           data=md, REML=F,
           control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))
summary(m1)

