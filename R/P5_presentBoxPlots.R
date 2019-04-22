setwd("E:/Dropbox/_PublishedWork/MatrixCompletion/myAirCoachDataCompletionFitbitSpire/R")

rm(list = ls())

absoluteError <- function(Ximp,Xreal){
  return(abs(Xreal-Ximp))
}

relativeError<- function(Ximp,Xreal){
  return(100*abs(Xreal-Ximp)/Xreal)
}


library(missForest)
library(mixtools)
library(oce)
library(ggplot2) 


datasetid<-96;
strid<-sprintf("%03d",datasetid)

datasetfolder<-"../_CSV/";
dataset<-paste("ProcessedDataset",strid,"/",sep=""); 
resultsfolder<-"../_CSV/Results/";
datasettarget<-paste("dataset",strid,sep="");  




yKNN<-data.frame()
fnm<-paste(resultsfolder,datasettarget,'_knn.csv',sep="")
M1=as.matrix(read.csv(fnm,header = FALSE))
#M1[,2]<-abs(M1[,2]-InitOutput)
#M1[,2]<-absoluteError(M1[,2],InitOutput)
#M1[,2]<-relativeError(M1[,2],InitOutput)
i<-5
for(i in seq(5,40,by=5)){
  rn<-paste('V',as.character(i),sep="")
  df<-data.frame(V= M1[which(M1[,1]==i),3] )
  names(df)[1]=rn
yKNN<-rbind(yKNN,t(df) )
}
yKNN<-t(yKNN)
yKNN<-as.data.frame(yKNN)








yLMC<-data.frame()
fnm<-paste(resultsfolder,datasettarget,'_mc_laplacian.csv',sep="")
M2=as.matrix(read.csv(fnm,header = FALSE))
#M2[,2]<-abs(M2[,2]-InitOutput)
#M2[,2]<-absoluteError(M2[,2],InitOutput)
#M2[,2]<-relativeError(M2[,2],InitOutput)
for(i in seq(5,40,by=5)){
  rn<-paste('V',as.character(i),sep="")
  df<-data.frame(V= M2[which(M2[,1]==i),3] )
  names(df)[1]=rn
  yLMC<-rbind(yLMC,t(df) )
}
yLMC<-t(yLMC)
yLMC<-as.data.frame(yLMC)









yMC<-data.frame()
fnm<-paste(resultsfolder,datasettarget,'_mc.csv',sep="")
M3=as.matrix(read.csv(fnm,header = FALSE))
#M3[,2]<-abs(M3[,2]-InitOutput)
#M3[,2]<-absoluteError(M3[,2],InitOutput)
#M3[,2]<-relativeError(M3[,2],InitOutput)
for(i in seq(5,40,by=5)){
  rn<-paste('V',as.character(i),sep="")
  df<-data.frame(V= M3[which(M3[,1]==i),3] )
  names(df)[1]=rn
  yMC<-rbind(yMC,t(df) )
}
yMC<-t(yMC)
yMC<-as.data.frame(yMC)









yMF<-data.frame()
fnm<-paste(resultsfolder,datasettarget,'_missForest.csv',sep="")
M4=as.matrix(read.csv(fnm,header = FALSE))
#M4[,2]<-abs(M4[,2]-InitOutput)
#M4[,2]<-absoluteError(M4[,2],InitOutput)
#M4[,2]<-relativeError(M4[,2],InitOutput)
for(i in seq(5,40,by=5)){
  rn<-paste('V',as.character(i),sep="")
  df<-data.frame(V= M4[which(M4[,1]==i),3] )
  names(df)[1]=rn
  yMF<-rbind(yMF,t(df) )
}
yMF<-t(yMF)
yMF<-as.data.frame(yMF)








M1<-as.data.frame(M1)
M1<-cbind(M1[,cbind(1,3)],'KNN')
colnames(M1)<-cbind('V1','V2','V3')
M2<-as.data.frame(M2)
M2<-cbind(M2[,cbind(1,3)],'Laplacian Matrix Completion')
colnames(M2)<-cbind('V1','V2','V3')
M3<-as.data.frame(M3)
M3<-cbind(M3[,cbind(1,3)],'Matrix Completion')
colnames(M3)<-cbind('V1','V2','V3')
M4<-as.data.frame(M4)
M4<-cbind(M4[,cbind(1,3)],'MissForest')
colnames(M4)<-cbind('V1','V2','V3')






MM<-as.data.frame(rbind(M1,M2,M3,M4))
df<-MM
MM$V1<-factor(MM$V1)







GRAPH<-ggplot(MM, aes(x = V1, y = V2, fill = V3)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "NRMSE",breaks = seq(0, 1, 0.1),limits=c(0, 0.4)) +
  scale_x_discrete(name = "Missing entries %") +
  #ggtitle("Boxplot of mean ozone by month") +
  theme_bw() +
  theme(plot.title = element_text(size = 14, family = "Tahoma", face = "bold"),
        text = element_text(size = 32, family = "Tahoma",lineheight = 38),
        # axis.title = element_text(face=""),
        axis.text.x=element_text(size = 14),
        legend.position = c(0.22, 0.85),
        legend.direction = "vertical") +
  scale_fill_brewer(palette = "Accent") +
  labs(fill = "")


GRAPH

