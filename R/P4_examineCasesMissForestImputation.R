##MISSFOREST
setwd("E:/Dropbox/_PublishedWork/MatrixCompletion/myAirCoachDataCompletionFitbitSpire/R")
rm(list = ls())
library(doParallel)
library(missForest)



#datasetid<-32;
#datasetid<-31;
#datasetid<-33;

datasetid<-96

strid<-sprintf("%03d",datasetid)

datasetfolder<-"../_CSV/";
dataset<-paste("ProcessedDataset",strid,"/",sep=""); 
resultsfolder<-"../_CSV/Results/";
datasettarget<-paste("dataset",strid,sep="");  


writecompleted<-TRUE;
completedfolder<-paste(datasetfolder,dataset,'Completed/',sep="")

path<-paste(datasetfolder,dataset,sep="");
fnms<-list.files(path)
fnms<-fnms[4:(length(fnms))]
Data=read.table(paste(path,"initial.csv",sep=""),header = FALSE, sep = ",")
a <- array(0,dim=c(0,5))
i<-fnms[4]

registerDoParallel(cores=4)


for (i in fnms){
  #print(fnms)
  split1<-strsplit(i, "\\.")[[1]]
  split2<-strsplit(split1[1], "_")[[1]]
  missing<-as.numeric(split2[3])
  permnumber<-as.numeric(split2[4])
  Case=read.table(paste(path,i,sep=""),header = FALSE, sep = ",")
  #Case[Case==0]=NA
  Imp<-missForest(Case, maxiter = 2, ntree = 20,parallelize = "forests")
  #Imp<-missForest(Case)
  Result<-Imp$ximp
  Xreal<-data.matrix(Data)
  Ximp <-data.matrix(Result)
  
  if(writecompleted){
    write.table(Ximp,file=paste(completedfolder,"completed_MF_",i,sep=""),row.names=FALSE,col.names=FALSE, sep=",")
  }
  
  
  dImp<-Ximp[2:(length(Ximp))]-Ximp[1:(length(Ximp)-1)]
  dReal<-Xreal[2:(length(Xreal))]-Xreal[1:(length(Xreal)-1)]
  

  
  dPE<-norm(as.matrix(dReal-dImp))/norm(as.matrix(dReal))
  dNMSE<-(mean((dReal-dImp)^2)/mean((dReal-mean(dReal))^2))
  PE<-norm(Xreal-Ximp)/norm(Xreal)
  NMSE<-(mean((Xreal-Ximp)^2)/mean((Xreal-mean(Xreal))^2))
  NRMSE<-sqrt(NMSE)
  h1<-hist(as.vector(as.matrix(Xreal)),xlim=c(40.0,180.0),breaks=140)
  h2<-hist(as.vector(as.matrix(Ximp)),xlim=c(40.0,180.0),breaks=140)
  
  hG=h1$counts
  hG=as.numeric(hG)
  hG[is.na(hG)]=0
  hR=h2$counts
  hR=as.numeric(hR)
  hR[is.na(hR)]=0
  dHist<-100*(mean((hG-hR)^2)/mean((hG-mean(hG))^2))
  #dHist<-100*norm(as.matrix(hG-hR))/norm(as.matrix(hG))
  PE<-100*PE

  
  
  print(paste("Missing : ", missing,"%","  NMSE : ",NMSE,"%","  NRMSE : ",NRMSE,"%","   PE : ",PE,"%","  dNRMSE : ",dNMSE,"%","   dPE : ",dPE,"%","   dHist : ",dHist,"%",sep="" ))
  a<-rbind(a,c(missing,NMSE,NRMSE,PE,dHist))
}

b<-a[complete.cases(a),]
targetfilepath<-paste(resultsfolder,datasettarget,"_missForest.csv",sep="")
write.table(b,file=targetfilepath,row.names=FALSE,col.names=FALSE, sep=",")




