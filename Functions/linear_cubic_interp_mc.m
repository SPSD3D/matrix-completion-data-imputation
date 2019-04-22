clc;
close all;
clear variables;
MInit=csvread('CSV/testing.csv');

R1c=[];

R2c=[];
repetitions=4;
for ccc=1:repetitions

R1=[];
for missingEntriesAsAPercentage=0.05:0.05:0.75;
MM=MInit;
r = rand(size(MM));
MM(r>(1-missingEntriesAsAPercentage))=0;
MissingValuesIndices=find(MM==0);
MissingValues=MInit(MissingValuesIndices);
realMissing=100*numel(MM(MM==0))/numel(MM);
y=reshape(MM,1,[]);
nz=find(y);
y=[y(nz(1)) y y(nz(end))]; %padding
x=linspace(1,length(y),length(y));
xq=1:length(y);
x(y==0)=[];
y(y==0)=[];
MC = interp1(x,y,xq);
MC=MC(2:end-1);
xq=1:length(MC);
GroundTruth=reshape(MInit,1,[]);
GroundTruth(isnan(MC))=[];
MC(isnan(MC))=[];
%NMSE=(1.0/numel(GroundTruth))*sum(((GroundTruth-MC).^2)./(GroundTruth.^2));
NMSE=(mean((GroundTruth-MC).^2)/mean((GroundTruth-mean(GroundTruth)).^2));
R1=[R1;[realMissing 100*NMSE]];
end
if numel(R1c)==0
    R1c=R1;
else
    R1c=R1c+R1;
end





R2=[];
for missingEntriesAsAPercentage=0.05:0.05:0.75;
MM=MInit;
r = rand(size(MM));
MM(r>(1-missingEntriesAsAPercentage))=0;
MissingValuesIndices=find(MM==0);
MissingValues=MInit(MissingValuesIndices);
realMissing=100*numel(MM(MM==0))/numel(MM);
y=reshape(MM,1,[]);
nz=find(y);
y=[y(nz(1)) y y(nz(end))]; %padding
x=linspace(1,length(y),length(y));
xq=1:length(y);
x(y==0)=[];
y(y==0)=[];
MC = interp1(x,y,xq,'Spline');
MC=MC(2:end-1);
xq=1:length(MC);
GroundTruth=reshape(MInit,1,[]);
GroundTruth(isnan(MC))=[];
MC(isnan(MC))=[];
NMSE=(mean((GroundTruth-MC).^2)/mean((GroundTruth-mean(GroundTruth)).^2));

R2=[R2;[realMissing 100*NMSE]];
end
if numel(R2c)==0
    R2c=R2;
else
    R2c=R2c+R2;
end

end

R1c=R1c/repetitions;
R2c=R2c/repetitions;
figure(300)
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 16 9])
plot(R1c(:,1),R1c(:,2),'-o','LineWidth',2,'MarkerSize',4);
hold on;
plot(R2c(:,1),R2c(:,2),'-.','LineWidth',2,'MarkerSize',4);
hold off;
xlabel({'Missing entries uniformly distributed(%) '})
ylabel({'Error(%)'})
title('Error vs Missing entries')
legend({'Linear interpolation','Spline'},'FontSize',16)
saveas(300,[num2str(300) '.jpg'])
