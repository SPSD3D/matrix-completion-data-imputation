##MISSFOREST
setwd("E:/Projects/myAirCoachDataCompletionFitbitSpire/R")
rm(list = ls())
library(missForest)
library(mixtools)

datasetid<-2;

strid<-sprintf("%03d",datasetid)

datasetfolder<-"../CSV/";
dataset<-paste("ProcessedDataset",strid,"/",sep=""); 
resultsfolder<-"../CSV/Results/";
datasettarget<-paste("dataset",strid,sep="");  


writecompleted<-TRUE;
completedfolder<-paste(datasetfolder,dataset,'Completed/',sep="")

path<-paste(datasetfolder,dataset,sep="");
fnms<-list.files(path)
fnms<-fnms[4:(length(fnms))]
Data=read.table(paste(path,"initial.csv",sep=""),header = FALSE, sep = ",")

d=as.numeric(as.vector(as.matrix(Data)))
f=spectrum(d)
f=f$freq


mixmdl = normalmixEM(d,3)
plot(mixmdl,which=2)
lines(density(d), lty=2, lwd=2)


x=seq(1, 200, by = 0.5)
s<-mixmdl$sigma[1]
mu<-mixmdl$mu[1]
y1=(1/sqrt(2*pi*s))*exp(-((x-mu)^2)/(2*s^2))
s<-mixmdl$sigma[2]
mu<-mixmdl$mu[2]
y2=(1/sqrt(2*pi*s))*exp(-((x-mu)^2)/(2*s^2))

ytot<-mixmdl$lambda[1]*y1+mixmdl$lambda[2]*y2

plot(x,mixmdl$lambda[1]*y1+mixmdl$lambda[2]*y2)
