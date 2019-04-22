if !isdefined(:definitions)
  include("simpleMatrixCompletion.jl")
  include("polynomialMatrixCompletion.jl")
  include("proposed.jl")
  include("fns.jl")
  include("..\lib\MixtureModels\MixtureModels.jl")
  using Gadfly, Plots, ColorSchemes,Dierckx,Interpolations,StatsBase,Base
  definitions=true;
end



test=true;
uniform=false;
block=false;
uniform2=false;
checker=false;
doLaplacian=false;
ttype="row"
gamma=0.2;
mRank=4;
mStep=0.2;
mIter=1000;
printToFiles=true
replaceOrig=true


#Dataset
datasetid=1;

str=lpad(datasetid,3,0)

datasetfolder="../CSV/"
dataset="ProcessedDataset" * str * "/"
resultsfolder="../CSV/Results/"
datasettarget="dataset" * str


writecompleted=true;
completedfolder=datasetfolder*dataset*"Completed/";


Data = readdlm(datasetfolder*dataset*"initial.csv", ',', Float64);
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

v,bins=myhist(vec(Data),65.0,180.0,20)
Plots.plot(v)


for i=1:size(ListOfDataFiles,1)
#for i=60:61
  fnm=ListOfDataFiles[i];
  sbstrs=split(fnm,".")
  sbstrs=split(sbstrs[1],"_")
  missing=parse(Int64,sbstrs[3])
  permnumber=parse(Int64,sbstrs[4])



  Case = readdlm(datasetfolder*dataset*fnm, ',', Float64);
  #Case=transpose(Case)

  Omega=zeros(size(Case))
  Omega[find(Case)]=1.0

  Result=mc_reconstruction_simple(Data,Omega,mRank,0.1,mStep,mIter,false);
  if replaceOrig
    #inds=find(transpose(Omega))
    inds=find((Omega))
    Result[inds]=Data[inds]
  end

  if writecompleted
    fw = open(completedfolder * "completed_mc_"*fnm, "w")
    writedlm(fw,Result,",")
    close(fw)
  end




  dData = Data[2:(end)]-Data[1:(end-1)]
  dResult = Result[2:(end)]-Result[1:(end-1)]
  #dData=reshape(dData, size(Data,1),size(Data,2))
  #dResult=reshape(dResult, size(Result,1),size(Result,2))
  #Plots.plot(dData)
  #Plots.plot(dResult)
  #Plots.plot(dResult-dData)
  PE=vecnorm(Result-Data)/vecnorm(Data)
  NMSE=(mean((Result-Data).^2)/mean((Data-mean(Data)).^2))
  NRMSE=sqrt(NMSE)
  dPE=vecnorm(dResult-dData)/vecnorm(dData)
  dNMSE=(mean((dResult-dData).^2)/mean((dData-mean(dData)).^2))
  hG,bins=myhist(vec(Data),65.0,180.0,20)
  hR,bins=myhist(vec(Result),65.0,180.0,20)

  dHist=100*vecnorm(hR-hG)/vecnorm(hG)

  PE=PE*100
  NRMSE=NRMSE*100
  NMSE=NMSE*100
  dNMSE=dNMSE*100
  dPE=dPE*100
  println("Missing : " * string(missing) * "   NMSE : " * string((NMSE)) * "  " * "NRMSE : " * string((NRMSE)) *"  PE : " * string((PE)) *"  dPE : " * string((dPE)) * "  dNMSE : " * string((dNMSE)) * "  Hist : " * string((dHist)))
  resultsCollection=vcat(resultsCollection,[missing NMSE NRMSE PE dHist])
end
resultsCollection=resultsCollection[1:size(resultsCollection,1) .!= 1,: ]

f = open(resultsfolder * datasettarget * "_mc.csv", "w")
writedlm(f,resultsCollection,",")
close(f)
