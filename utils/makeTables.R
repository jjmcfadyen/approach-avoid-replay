makeTables <- function(m,mnum) {
# m is the model object
# mnum is the number of this model
  
  require(sjPlot)
  S <- summary(m)
  
  thistable <- as.data.frame(S$coefficients)
  thistable <- round(thistable,3)
  
  pvals <- select(thistable,starts_with("P"))
  thistable$star <- rep("",nrow(thistable))
  thistable$star[pvals<=0.05] = "*"
  thistable$star[pvals<=0.01] = "**"
  thistable$star[pvals<=0.001] = "***"
  if (any(names(thistable)=="df")){
    thistable <- thistable[,c(1,2,5,6)]
  } else {
    thistable <- thistable[,c(1,2,4,5)] # includes additional df variable
  }
  
  thistable$P[thistable$P=="0"] = "< 0.001"
  
  N1 <- c("RewVal","LossVal","RewProb","Replay_typerewarding","Choice1","PathRecency","PathProb","Replay_type","Replay_differential","EV","Replay_rewarding","Replay_aversive","PCA1","PCA2")
  N2 <- c("Reward magnitude","Loss magnitude","Transition probability","Replay type","Choice","Recency","Transition probability","Replay type","Differential replay","Expected value","Rewarding path replay","Punishing path replay","Risk-aversion","Anxiety")
  
  R <- row.names(thistable)
  for (i in 1:length(N1)){
    R <- gsub(N1[i],N2[i],R)
  }
  R <- gsub(":"," × ",R)
  thistable$Effect <- R
  thistable <- thistable[,c(6,1,2,5,4)]
  
  names(thistable)[1] <- "Fixed Effect"
  names(thistable)[2] <- "??"
  names(thistable)[3] <- "SEM"
  names(thistable)[4] <- "p"
  
  fname <- paste("modeltable_",mnum,".doc",sep="")
  
  mname <- as.character(formula(m))
  mname <- paste("Model ", mnum, ": ",mname[2]," ~ ",mname[3],sep="",collapse="")
  for (i in 1:length(N1)){
    mname <- gsub(N1[i],N2[i],mname)
  }
  mname <- gsub("\\*","×",mname)
  
  v <- round(vif(m),2)
  v <- paste("Variable inflation factor = ",min(v)," to ",max(v),sep="")
  
  autocorr <- durbinWatsonTest(resid(m))
  autocorr <- paste("Durbin-Watson = ",round(autocorr,2),sep="")
  
  tab_df(thistable,
         file=fname,
         title=mname,
         footnote=paste(v,autocorr,sep=", "),
         show.footnote=TRUE,
         show.rownames=FALSE)

}