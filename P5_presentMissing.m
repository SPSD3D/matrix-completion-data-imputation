clc;
clearvars -except id
close all;

if ~exist('id','var')
id=2;
end

idstr = sprintf('%03d',id);
missingprev=0;
datasetfolder='CSV/';
dataset=['ProcessedDataset' idstr '/'];
resultsfolder='CSV/Results/';
datasettarget=['dataset' idstr];
writecompleted=true;
completedfolder=[datasetfolder dataset 'Completed/'];
Data=csvread([datasetfolder dataset 'initial.csv']);
Omega=zeros(size(Data));
resultsCollection=[];
H=dir([datasetfolder dataset '*.csv']);
H=H(2:end,:);

ts=Data(:);
ts=ts';

imagesc(Data)

figure(10)
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 16 8])
ci=1;

for k=1:length(H)
%for k=10:20
   fnm=H(k).name;
   spl = strsplit(fnm,'.');
   spl2=strsplit(spl{1},'_');
   missing=str2num(spl2{3});
   permnumber=str2num(spl2{4});
   Case=csvread([datasetfolder dataset fnm]);
   y=reshape(Case,1,[]);
   nz=find(y);
   y=[y(nz(1)) y y(nz(end))]; %padding
   x=linspace(1,length(y),length(y));
   xq=1:length(y);
   x(y==0)=[];
   y(y==0)=[];
   Result = pchip(x,y,xq);
   Result=Result(2:end-1);
   xq=1:length(Result);
   GroundTruth=reshape(Data,1,[]);
   
 
   if missingprev~=missing
   subplot(2,4,ci)
   colormap(gray)
   imagesc(im2bw(Case, 0.5));
   axis off
   missingprev=missing;
   ci=ci+1;
   end
   
   
   if writecompleted
   completed=reshape(Result,size(Data));
   csvwrite([completedfolder  'completed_hermite_' fnm],completed);
   end
   
   %GroundTruth(isnan(Result))=[];
   %Result(isnan(Result))=[];
   
   deltaGroundTruth=GroundTruth(2:(end))-GroundTruth(1:(end-1));
   deltaResult=Result(2:(end))-Result(1:(end-1));
   PE=norm(GroundTruth-Result)/norm(GroundTruth);
   NMSE=(mean((GroundTruth-Result).^2)/mean((GroundTruth-mean(GroundTruth)).^2));
   NRMSE=sqrt(NMSE);
   PEDelta=norm(deltaGroundTruth-deltaResult)/norm(deltaGroundTruth);
   NMSEDelta=(mean((deltaGroundTruth-deltaResult).^2)/mean((deltaGroundTruth-mean(deltaGroundTruth)).^2));
   dHist=0;
%    Bins=linspace(40.5,179.5,140);
%    Ghist=histogram(GroundTruth,Bins);
%    GhistV=Ghist.Values;
%    GhistV=GhistV/sum(GhistV);
%    Rhist=histogram(Result,Bins);
%    RhistV=Rhist.Values;
%    RhistV=RhistV/sum(RhistV);
%    dHist=100*(mean((GhistV-RhistV).^2)/mean((GhistV-mean(GhistV)).^2));
      
   %GFreq=10*log10(abs(fft(GroundTruth).^2));
   %GFreq=GFreq(1:((length(GFreq)/2)+1));
   %RFreq=10*log10(abs(fft(Result).^2));
   %RFreq=RFreq(1:((length(RFreq)/2)+1));

%    figure(3)
%    plot(GhistV)
%    hold on
%    plot(RhistV)
%    hold off
   
   PE=100*PE;
   PEDelta=100*PEDelta;
   NMSE=100*NMSE;
   NMSEDelta=100*NMSEDelta;
   NRMSE=100*NRMSE;
   resultsCollection=[resultsCollection; [missing NMSE NRMSE PE dHist] ];
   disp(['Missing : ' num2str(missing)  '    NMSE : ' num2str(NMSE) '    NRMSE : ' num2str(NRMSE) '    PE : ' num2str(PE) '   dPE : ' num2str(PEDelta)  '   dNMSE : ' num2str(NMSEDelta) '   Hist : ' num2str(dHist) ]);
end
%csvwrite([resultsfolder  datasettarget '_hermite.csv'],resultsCollection);