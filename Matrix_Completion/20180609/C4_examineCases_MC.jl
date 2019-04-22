if !isdefined(:definitions)
  include("definitions.jl")
  using Gadfly, Plots, ColorSchemes,Dierckx,Interpolations,StatsBase,Base
  definitions=true;
end
include("proposed.jl")

#Dataset & Path
datasetid=96;
str=lpad(datasetid,3,0)
datasetfolder="../_CSV/"
dataset="ProcessedDataset" * str * "/"
resultsfolder="../_CSV/Results/"
datasettarget="dataset" * str

mRank=4;
mStep=0.5;
mIter=400;
printToFiles=true
replaceOrig=true






writecompleted=true;
completedfolder=datasetfolder*dataset*"Completed/";
DataIn = readdlm(datasetfolder*dataset*"initial.csv", ',', Float64);
Udata,Sdata,Vdata=svd(DataIn)

#q,dq=make_quantizer(DataIn,16)
#M=float(q(DataIn))
#L,S=rpca(DataIn)


GroundTruth=deepcopy(DataIn)
Data=deepcopy(DataIn)




Case=zeros(size(Data))
Omega=zeros(size(Data))
Result=zeros(size(Data))
dData=zeros(size(Data))
dResult=zeros(size(Data))
ListOfDataFiles=readdir(datasetfolder*dataset)
deleteat!(ListOfDataFiles, findin(ListOfDataFiles, ["initial.csv"]))
deleteat!(ListOfDataFiles, findin(ListOfDataFiles, ["Completed"]))
deleteat!(ListOfDataFiles, findin(ListOfDataFiles, ["_readme.txt"]))
resultsCollection=[0.0 0.0 0.0 0.0 0.0];
for i=1:size(ListOfDataFiles,1)
#for i=1:3
  fnm=ListOfDataFiles[i];
  sbstrs=split(fnm,".")
  sbstrs=split(sbstrs[1],"_")
  missing=parse(Int64,sbstrs[3])
  permnumber=parse(Int64,sbstrs[4])
  Case = readdlm(datasetfolder*dataset*fnm, ',', Any);


  indices = find( x->(x != "NA"), Case)
  Omega=zeros(size(Case))
  Omega[indices]=1.0
  #Omega[find(Case)]=1.0



  Result,conv,iter=mc_reconstruction_simple(Data,Omega,mRank,0.1,mStep,mIter,false);

  if replaceOrig
    inds=find((Omega))
    Result[inds]=GroundTruth[inds]
  end
  if writecompleted
    fw = open(completedfolder * "completed_mc_"*fnm, "w")
    writedlm(fw,Result,",")
    close(fw)
  end
  PE=vecnorm(Result-GroundTruth)/vecnorm(GroundTruth)
  NMSE=(mean((Result-GroundTruth).^2)/mean((GroundTruth-mean(GroundTruth)).^2))
  NRMSE=sqrt(NMSE)


  hG,binsG=myhist(vec(GroundTruth),40.0,180.0,140)
  hG=hG/sum(hG)
  Plots.plot(hG)
  hR,binsR=myhist(vec(Result),40.0,180.0,140)
  hR=hR/sum(hR)
  Plots.plot(hR)
  dHist=100*(mean((hR-hG).^2)/mean((hG-mean(hG)).^2))


  dData = GroundTruth[2:(end)]-GroundTruth[1:(end-1)]
  dResult = Result[2:(end)]-Result[1:(end-1)]
  dPE=vecnorm(dResult-dData)/vecnorm(dData)
  dNMSE=(mean((dResult-dData).^2)/mean((dData-mean(dData)).^2))
  PE=PE*100

  println("Missing : " * string(missing) * "   NMSE : " * string((NMSE)) * "  " * "NRMSE : " * string((NRMSE)) *"  PE : " * string((PE)) *"  dPE : " * string((dPE)) * "  dNMSE : " * string((dNMSE)) * "  Hist : " * string((dHist)))
  resultsCollection=vcat(resultsCollection,[missing NMSE NRMSE PE dHist])
end
resultsCollection=resultsCollection[1:size(resultsCollection,1) .!= 1,: ]

f = open(resultsfolder * datasettarget * "_mc.csv", "w")
writedlm(f,resultsCollection,",")
close(f)
