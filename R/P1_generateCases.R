setwd("E:/Dropbox/_GroundWork/Release/myAirCoachDataCompletionFitbitSpire/R")
rm(list = ls())
library(missForest)

#
# Select dataset
#
dataset<-"ProcessedDataset051";
readfrom<-"_single_day_1_per_1_minute_in.csv";
permnum<-30
mstart<-0.05
mend<-0.40
mstep<-0.05

args = commandArgs(trailingOnly=TRUE)


if (length(args)==6) {
  dataset<-args[1];
  readfrom<-args[2];
  permnum<-as.numeric(args[3]);
  mstart<-as.numeric(args[4]);
  mend<-as.numeric(args[5]);
  mstep<-as.numeric(args[6]);
}









data=read.table(paste("../Data/",readfrom,sep=""),header = FALSE, sep = ",")
s <- svd(data)
D <- diag(s$d)
plot(s$d, type = "s", main = "Principal components")






path<-paste("../CSV/",dataset,"/initial.csv",sep="");
write.table(data,file=path,row.names=FALSE,col.names=FALSE, sep=",")






for (perm in 1:permnum){
  for (i in seq(mstart, mend, by = mstep)){
    data.mis <- prodNA(data, noNA = i)
    data.miszeros<- replace(data.mis, is.na(data.mis), 0)
    data.miszeros[1,1]=data[1,1]
    data.miszeros[nrow(data),ncol(data)]=data[nrow(data),ncol(data)]
    realmissing<-sum(is.na(data.mis))/(ncol(data.mis)*nrow(data.mis))
    missingFilename<-paste("uniform_missing",sprintf("%02d",trunc(100*i)),sprintf("%04d", perm),sep="_",collapse = NULL)
    write.table(data.miszeros,file=paste("../CSV/",dataset,"/",missingFilename,".csv",sep=""),row.names=FALSE,col.names=FALSE, sep=",")
  }
}


