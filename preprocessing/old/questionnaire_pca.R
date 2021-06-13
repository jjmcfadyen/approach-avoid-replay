library(dplyr)
library(factoextra)

setwd("D:/Files/2020_RiskyReplay/2021")

Q <- read.csv("groupBehaviour.csv")

qpca <- Q %>% 
  select(-Subject) %>%
  prcomp(center = TRUE,scale. = TRUE)
qpca.summary <- summary(qpca)

fviz_eig(qpca) # scree plot
fviz_pca_var(qpca, # dimension loadings
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#8700FF", "#D800FF", "#F74F8C", "#FF851E", "#FFD800"),
             repel = TRUE     # Avoid text overlapping
)

# pick components with > 1 eigenvalue
qpca.comp <- qpca$x[,qpca.summary$importance[1,] > 1] %>% as_tibble() %>% mutate(Subject = unique(Q$Subject))

# add the PCA components to the main dataframe
for (p in 1:(dim(qpca.comp)[2]-1)) {
  vname <- paste("PC",as.character(p),sep="")
  Q[[vname]] = rep(NA,dim(Q)[1])
  for (subject in unique(Q$Subject)) {
    Q[[vname]][Q$Subject==subject] = qpca.comp[[vname]][qpca.comp$Subject==subject]
  }
}

# plot mappings
longcomp <- gather(Q, PCA, PScore, PC1:PC2, factor_key=T)
longcomp <- gather(longcomp,Questionnaire,QScore,IUS:Social)

ggplot(data=longcomp, aes(x=QScore,y=PScore,color=PCA)) + 
  geom_point() + 
  geom_smooth(method='lm', se=FALSE) + 
  facet_wrap(~PCA+Questionnaire) + 
  theme_classic() 

# save data
Qsave <- Q
names(Qsave)[9:10] <- c("RiskSeeking","Anxiety")
Qsave$Anxiety <- Qsave$Anxiety * (-1)

write.csv("questionnaire_PCA.csv",row.names=F)
