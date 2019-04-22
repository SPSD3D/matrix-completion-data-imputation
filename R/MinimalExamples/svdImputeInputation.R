##impute.svd
setwd("E:/Projects/myAirCoachDataCompletionFitbitSpire/R")

rm(list = ls())
library(bcv)
data=read.csv("../CSV/testing.csv", fill = TRUE)
a <- array(0,dim=c(0,2))
for (i in seq(0.05, 0.35, by = 0.05)){
data.mis <- prodNA(data, noNA = i)
data.imp <-impute.svd(data.mis, k = 2,  maxiter=10000)$x
data.impf<-data.imp
Xreal<-data.matrix(data)
Ximp <-data.matrix(data.impf)
err<-(mean((Xreal-Ximp)^2)/mean((Xreal-mean(Xreal))^2))

#print(paste("Missing:", i,"Error:",data.imp$OOBerror))

print(paste("Missing:", 100*i,"%","  ","Error:",100*err,"%"))
a<-rbind(a,c(100*i,100*err))
}
plot(a[,1],a[,2],type = "o")
