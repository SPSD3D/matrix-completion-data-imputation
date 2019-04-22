cd("E:/Dropbox/_PublishedWork/MatrixCompletion/myAirCoachDataCompletionFitbitSpire/Matrix_Completion/20180609")

#if !isdefined(:definitions)
  include("definitions.jl")
  using Gadfly, Plots, ColorSchemes,Dierckx,Interpolations,StatsBase,Base
  definitions=true;
#end
include("proposed.jl")


#Dataset & Path
datasetid=96;
gamma=0.4;
mRank=20;
mStep=0.4;
mIter=600;

str=lpad(datasetid,3,0)
datasetfolder="../../_CSV/"
dataset="ProcessedDataset" * str * "/"
resultsfolder="../../_CSV/Results/"
datasettarget="dataset" * str

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
resultsCollection=[0.0 0.0 0.0 0.0 0.0 0.0 0.0];
#for i=1:size(ListOfDataFiles,1)
A=spzeros(size(Data,1),size(Data,1));

convergence_curve_mc=0;
convergence_curve_mc_sum=0;
error_between_iterations=0;
error_between_iterations_sum=0;
iter=0;

#Lc=dimensionLaplacian(Data,1,1)
#Lr=dimensionLaplacian(Data,1,1)
#L=interpL(vec(Data))

  directorylen=size(ListOfDataFiles,1)

  for i=1:directorylen
  #for i=1:5
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
  indsz=find( x->(x == "NA"),Case)
  L1=interpL(Case,indsz)
  #L=interpL(vec(vec(Data).*vec(Omega)))
  #L=sparse(L)



  #F=fft(Data[:])
  #F=abs(F).^2
  #Plots.plot(F[2:trunc(Int32,(length(F)/2)+1)])
  #L2=interpL(vec(vec(transpose(Data)).*vec(transpose(Omega))))
  #L2=sparse(L2)
  #Result,convergence_curve_mc,iter,error_between_iterations=mc_reconstruction_laplacian_interpolation(Data,Omega,L1,L2,mRank,0.1,mStep,mIter,false);



  tic()
  Result,convergence_curve_mc,iter,error_between_iterations=mc_reconstruction_laplacian_interpolation(Data,Omega,L1,mRank,0.1,mStep,mIter,false,gamma);
  timed=toc()

  if i==1
    convergence_curve_mc_sum=convergence_curve_mc
    error_between_iterations_sum=error_between_iterations
  else
    convergence_curve_mc_sum=convergence_curve_mc_sum+convergence_curve_mc
    error_between_iterations_sum=error_between_iterations_sum+error_between_iterations
  end


  if replaceOrig
    inds=find((Omega))
    Result[inds]=GroundTruth[inds]
  end
  if writecompleted
    fw = open(completedfolder * "completed_LMC_"*fnm, "w")
    writedlm(fw,Result,",")
    close(fw)
  end
  PE=vecnorm(Result-GroundTruth)/vecnorm(GroundTruth)
  NMSE=(mean((Result-GroundTruth).^2)/mean((GroundTruth-mean(GroundTruth)).^2))
  NRMSE=sqrt(NMSE)


  #fG=10log10(abs(fft(vec(GroundTruth))).^2)
  #fG=fG[1:trunc(Int,(length(fG)/2)+1)]
  #Plots.plot(fG)
  #fR=10log10(abs(fft(vec(Result))).^2)
  #fR=fR[1:trunc(Int,(length(fR)/2)+1)]
  #Plots.plot(fR)
  #U,Sd,V=svd(reshape(fG,size(GroundTruth)))


  hG,binsG=myhist(vec(GroundTruth),40.0,180.0,140)
  hG=hG/sum(hG)
  Plots.plot(hG)
  hR,binsR=myhist(vec(Result),40.0,180.0,140)
  hR=hR/sum(hR)
  Plots.plot(hR)
  dHist=100*(mean((hR-hG).^2)/mean((hG-mean(hG)).^2))

  #dHist=100*vecnorm(hR-hG)/vecnorm(hG)

  dData = GroundTruth[2:(end)]-GroundTruth[1:(end-1)]
  dResult = Result[2:(end)]-Result[1:(end-1)]
  dPE=vecnorm(dResult-dData)/vecnorm(dData)
  dNMSE=(mean((dResult-dData).^2)/mean((dData-mean(dData)).^2))


  PE=PE*100
  #NRMSE=NRMSE*100
  #NMSE=NMSE*100
  #dNMSE=dNMSE*100
  dPE=dPE*100
  println("Missing : " * string(missing) * "   NMSE : " * string((NMSE)) * "  " * "NRMSE : " * string((NRMSE)) *"  PE : " * string((PE)) *"  dPE : " * string((dPE)) * "  dNMSE : " * string((dNMSE)) * "  Hist : " * string((dHist)) * "  # of iter : " * string((iter)) )
  resultsCollection=vcat(resultsCollection,[missing NMSE NRMSE PE dHist iter timed])
end
resultsCollection=resultsCollection[1:size(resultsCollection,1) .!= 1,: ]

f = open(resultsfolder * datasettarget * "_mc_laplacian.csv", "w")
writedlm(f,resultsCollection,",")
close(f)



convergence_curve_mc_sum=convergence_curve_mc_sum/directorylen
error_between_iterations_sum=error_between_iterations_sum/directorylen



Plots.plot(convergence_curve_mc_sum[1:30],label=["NRMSE vs number of iteration" ""],lw=3)
xlabel!("Iteration")
ylabel!("NRMSE")



Plots.plot(error_between_iterations_sum[1:30],label=["Err vs number of iteration" ""],lw=3)
xlabel!("Iteration")
ylabel!("Err")


Plots.plot(resultsCollection[2:100,7],label=["NRMSE vs number of iteration" ""],lw=3)
xlabel!("Iteration")
ylabel!("NRMSE")
