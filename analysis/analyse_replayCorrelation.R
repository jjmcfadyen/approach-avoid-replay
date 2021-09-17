library(lme4)
library(lmerTest)
library(interactions)
library(car)
library(emmeans)
library(tidyverse)
library(viridis)
library(extrafont)
library(gridExtra)

setwd("D:/2020_RiskyReplay/results")
d <- read.csv("replay_correlation.csv")

# Add info
d$Catch <- (d$nV_1>0 & d$nV_2>0) | (d$nV_1<0 & d$nV_2<0)

d$P <- round(d$P,1)
d$Certainty <- d$P
d$Certainty[d$P==.1] <- .9
d$Certainty[d$P==.3] <- .7

d$BinChoice <- d$Choice
d$Choice[d$Choice==1] <- "approach"
d$Choice[d$Choice==2] <- "avoid"
d$Choice <- as.factor(d$Choice)

# QUESTIONNAIRES
Q <- read.csv("D:/2020_RiskyReplay/results/questionnaire_results.csv")
Q <- Q[Q$Subject %in% intersect(unique(d$Subject),Q$Subject),]

pca <- prcomp(Q[,c(2,3,5,6,7,8,9)],center=TRUE,scale.=TRUE)
eig <- summary(pca)
PC <- as.data.frame(pca$x)

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

# LME
md <- d %>%
  filter(Forced==0,Catch==FALSE,PathIteration==1)
md$RT <- scale(md$RT,center=TRUE,scale=FALSE)
md$Certainty <- scale(md$Certainty,center=TRUE,scale=FALSE)

m0 <- lmer(CrossLag ~ Choice*Certainty + (1|Subject),
            data=md, REML=FALSE,
            control=lmerControl(optimizer="bobyqa",optCtrl=list(maxfun=10e6)))
summary(m0)

pred <- jtools::make_predictions(m0,pred="Certainty",
                                 at=list(Choice=unique(md$Choice)),
                                 pred.values=unique(md$Certainty))

ss <- sim_slopes()
