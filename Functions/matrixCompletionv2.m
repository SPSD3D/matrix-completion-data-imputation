clc;
close all;
clear variables;
MInit=csvread('CSV/testing.csv');

R1=[];
for missingEntriesAsAPercentage=10:10:80;
MM=MInit;
r = randi([1,100],size(MM));
MM(r>100-missingEntriesAsAPercentage)=0;
MissingValuesIndices=find(MM==0);
MissingValues=MInit(MissingValuesIndices);
realMissing=100*numel(MM(MM==0))/numel(MM);

y=reshape(MM,1,[]);
nz=find(y)
y=[y(nz(1)) y y(nz(end))]; %padding
x=linspace(1,length(y),length(y));
xq=1:length(y);
x(y==0)=[];
y(y==0)=[];
MC = interp1(x,y,xq);
MC=MC(2:end-1)
xq=1:length(MC);

GroundTruth=reshape(MInit,1,[]);
figureID=missingEntriesAsAPercentage+100;
figure(figureID)
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 16 9])
subplot(3,4,[1:3 5:7 9:11] )
plot(xq,GroundTruth,'--b','LineWidth',3)
hold on;
plot(xq,MC,'-og','LineWidth',1,'MarkerSize',2)
hold on;
plot(MissingValuesIndices,MissingValues,'xr','LineWidth',1,'MarkerSize',7)
hold off;
legend({'Ground truth','Matrix completion'},'FontSize',14)
subplot(3,4,4 )
imagesc(MM);
GroundTruth(isnan(MC))=[];
MC(isnan(MC))=[];
NMSE=(1.0/numel(GroundTruth))*sum(((GroundTruth-MC).^2)./(GroundTruth.^2));



subplot(3,4,12)
txt1 = ['Missing Entries:' num2str(realMissing) '%'];
txt2 = ['NMSE:' num2str(100*NMSE) '%'];
axis off
text(0.2,0.2,txt2,'FontSize',20)
text(0.2,0.7,txt1,'FontSize',20)
saveas(figureID,[num2str(figureID) '.jpg'])
R1=[R1;[realMissing 100*NMSE]];
end





R2=[];
for missingEntriesAsAPercentage=10:10:80;
MM=MInit;
r = randi([1,100],size(MM));
MM(r>100-missingEntriesAsAPercentage)=0;
MissingValuesIndices=find(MM==0);
MissingValues=MInit(MissingValuesIndices);
realMissing=100*numel(MM(MM==0))/numel(MM);

y=reshape(MM,1,[]);
nz=find(y)
y=[y(nz(1)) y y(nz(end))]; %padding
x=linspace(1,length(y),length(y));
xq=1:length(y);
x(y==0)=[];
y(y==0)=[];
MC = interp1(x,y,xq,'Spline');
MC=MC(2:end-1)
xq=1:length(MC);

GroundTruth=reshape(MInit,1,[]);
figureID=missingEntriesAsAPercentage+200;
figure(figureID);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 16 9])
subplot(3,4,[1:3 5:7 9:11] )
plot(xq,GroundTruth,'--b','LineWidth',3)
hold on;
plot(xq,MC,'-og','LineWidth',1,'MarkerSize',2)
hold on;
plot(MissingValuesIndices,MissingValues,'xr','LineWidth',1,'MarkerSize',7)
hold off;
legend({'Ground truth','Matrix completion'},'FontSize',14)
subplot(3,4,4)
imagesc(MM);
GroundTruth(isnan(MC))=[];
MC(isnan(MC))=[];
NMSE=(1.0/numel(GroundTruth))*sum(((GroundTruth-MC).^2)./(GroundTruth.^2));
%NMSE=norm(GroundTruth-MC,2)/norm(GroundTruth,2);
subplot(3,4,12)
txt1 = ['Missing Entries:' num2str(realMissing) '%'];
txt2 = ['NMSE:' num2str(100*NMSE) '%'];
axis off
text(0.2,0.2,txt2,'FontSize',20)
text(0.2,0.7,txt1,'FontSize',20)
saveas(figureID,[num2str(figureID) '.jpg'])
R2=[R2;[realMissing 100*NMSE]];
%sum(abs(GroundTruth-MC)./abs(GroundTruth))
end




R0=csvread('CSV/nmse_vs_missing_uniform.csv');
figure(300)
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 16 9])
plot(R0(:,1),R0(:,2),'-o','LineWidth',2,'MarkerSize',4);
hold on;
plot(R1(:,1),R1(:,2),'-o','LineWidth',2,'MarkerSize',4);
hold on;
plot(R2(:,1),R2(:,2),'-o','LineWidth',2,'MarkerSize',4);
hold off;
ylim([0 10]);
xlabel({'Missing entries uniformly distributed(%) '})
ylabel({'NMSE(%)'})
title('NMSE vs Missing entries for a single day heart rate measurments')
legend({'Nuclear norm minimization','Linear interpolation','Spline'},'FontSize',16)
saveas(300,[num2str(300) '.jpg'])
