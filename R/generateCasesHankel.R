setwd("E:/Dropbox/_GroundWork/Release/myAirCoachDataCompletionFitbitSpire/R")
rm(list = ls())
library(missForest)
source("functions.R")

permnum<-30
mstart<-0.05
mend<-0.40
mstep<-0.05

#Average
ma <- function(x,n=5){filter(x,rep(1/n,n), sides=2)}
#Median
mmed <- function(x,n=5){runmed(x,n)} 

readfrom<-"_single_day_1_per_1_minute_in.csv";

dataset<-"ProcessedDataset050";
data=read.table(paste("../Data/",readfrom,sep=""),header = FALSE, sep = ",")
s <- svd(data)
D <- diag(s$d)
plot(s$d, type = "s", main = "Principal components")









data2<-as.vector(matrix(as.matrix(data),nrow=1,byrow = FALSE))
data<-data2[1:128]
hankel<-as.matrix(matrix(as.matrix(data),nrow=1,byrow = FALSE))
b<-hnk(as.matrix(hankel),length(hankel)/2)
path<-paste("../CSV/",dataset,"/initial.csv",sep="");
write.table(b,file=path,row.names=FALSE,col.names=FALSE, sep=",")


data<-as.matrix(matrix(as.matrix(data),nrow=8,byrow = FALSE))

for (perm in 1:permnum){
  for (i in seq(mstart, mend, by = mstep)){
    data.mis <- prodNA(data, noNA = i)
    data.miszeros<- replace(data.mis, is.na(data.mis), 0)
    data.miszeros[1,1]=data[1,1]
    data.miszeros[nrow(data),ncol(data)]=data[nrow(data),ncol(data)]
    
    realmissing<-sum(is.na(data.mis))/(ncol(data.mis)*nrow(data.mis))
    
    hankel<-as.matrix(matrix(as.matrix(data.miszeros),nrow=1,byrow = FALSE))
    b<-hnk(as.matrix(hankel),length(hankel)/2)
    
    
    missingFilename<-paste("uniform_missing",sprintf("%02d",trunc(100*i)),sprintf("%04d", perm),sep="_",collapse = NULL)
    write.table(b,file=paste("../CSV/",dataset,"/",missingFilename,".csv",sep=""),row.names=FALSE,col.names=FALSE, sep=",")
  }
}


