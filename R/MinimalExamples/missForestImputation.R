##MISSFOREST
setwd("E:/Projects/myAirCoachDataCompletionFitbitSpire/R")
rm(list = ls())
library(missForest)
data=read.csv("../CSV/testing.csv", fill = TRUE)
counter<-1;
a <- array(0,dim=c(0,2))
for (i in seq(0.05, 0.75, by = 0.05)){
data.mis <- prodNA(data, noNA = i)
realmissing<-sum(is.na(data.mis))/(ncol(data.mis)*nrow(data.mis))
data.imp <- missForest(data.mis)
data.impf<-data.imp$ximp
Xreal<-data.matrix(data)
Ximp <-data.matrix(data.impf)
err<-(mean((Xreal-Ximp)^2)/mean((Xreal-mean(Xreal))^2))
#print(paste("Missing:", i,"Error:",data.imp$OOBerror))
print(paste("Missing:", 100*realmissing,"%","  ","Error:",100*err,"%"))
a<-rbind(a,c(100*realmissing,100*err))
}
plot(a[,1],a[,2],type = "o")
