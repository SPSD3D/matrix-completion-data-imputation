clc;
clear;
close all;


id=96;
idstr = sprintf('%03d',id);
resultToPresent=3
dataset=['dataset' idstr];

fileMC=['_CSV/Results/' dataset  '_mc.csv'];
fileLN=['_CSV/Results/' dataset '_linear.csv']
fileSPL=['_CSV/Results/' dataset '_spline.csv']
fileH=['_CSV/Results/' dataset '_hermite.csv']
fileMF=['_CSV/Results/' dataset '_missForest.csv']
fileGMC=['_CSV/Results/' dataset '_mc_laplacian.csv']
fileKNN=['_CSV/Results/' dataset '_knn.csv']
doMC=false;
doLN=false;
doSPL=false;
doH=false;
doMF=false;
doGMC=false;
doKNN=false;


part=5:5:50

if exist(fileKNN, 'file') == 2
Result_KNN=csvread(fileKNN);
Result_KNN(:,2:3)=Result_KNN(:,2:3)/100;
Result_KNN_mean=[];
doKNN=true;
for i=part
    Result_KNN_mean=[Result_KNN_mean;mean(Result_KNN(Result_KNN(:,1)==i,:),1)];
end
end

if exist(fileMC, 'file') == 2
Result_MC=csvread(fileMC);
Result_MC(:,2:3)=Result_MC(:,2:3)/100;
Result_MC_mean=[];
doMC=true;
for i=part
    Result_MC_mean=[Result_MC_mean;mean(Result_MC(Result_MC(:,1)==i,:),1)];
end
end

if exist(fileLN, 'file') == 2
Result_Linear=csvread(fileLN);
Result_Linear(:,2:3)=Result_Linear(:,2:3)/100;
Result_Linear_mean=[];
doLN=true;
for i=part
    Result_Linear_mean=[Result_Linear_mean;mean(Result_Linear(Result_Linear(:,1)==i,:),1)];
end
end

if exist(fileSPL, 'file') == 2
Result_Spline=csvread(fileSPL);
Result_Spline(:,2:3)=Result_Spline(:,2:3)/100;
Result_Spline_mean=[];
doSPL=true;
for i=part
    Result_Spline_mean=[Result_Spline_mean;mean(Result_Spline(Result_Spline(:,1)==i,:),1)];
end
end

if exist(fileH, 'file') == 2
Result_Hermite=csvread(fileH);
Result_Hermite(:,2:3)=Result_Hermite(:,2:3)/100;
Result_Hermite_mean=[];
doH=true;
for i=5:5:60
    Result_Hermite_mean=[Result_Hermite_mean;mean(Result_Hermite(Result_Hermite(:,1)==i,:),1)];
end
end

if exist(fileMF, 'file') == 2
Result_missForest=csvread(fileMF);
Result_missForest(:,2:3)=Result_missForest(:,2:3)/100;
Result_missForest_mean=[];
doMF=true;
for i=part
    Result_missForest_mean=[Result_missForest_mean;mean(Result_missForest(Result_missForest(:,1)==i,:),1)];
end
end

if exist(fileGMC, 'file') == 2
Result_MC_L=csvread(fileGMC);
Result_MC_L(:,2:3)=Result_MC_L(:,2:3)/100;
Result_MC_L_mean=[];
doGMC=true;
for i=part
    Result_MC_L_mean=[Result_MC_L_mean;mean(Result_MC_L(Result_MC_L(:,1)==i,:),1)];
end
executionTime_MC_L=Result_MC_L(:,6);
figure(20)
plot(executionTime_MC_L,'-o')
end


doMC=true;
doLN=false;
doSPL=false;
doH=false;
doMF=true;
doGMC=true;
doKNN=true;

showPE=false;
showHist=false;
limY.active1=false;
limY.active2=false;
limY.active3=false;
limY.value2=0.005;


line1.style='k-';
line1.width=2;

line2.style='k-s';
line2.width=2;

line3.style='k-x';
line3.width=2;

line4.style='k-*';
line4.width=2;

figure(1)
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 8.3 11.7])
if doMC
subplot(5,2,1)
plot(Result_MC(:,1),Result_MC(:,2),'b*');
title('Matrix Completion')
if limY.active1
ylim([0,100])
end
ylabel('NMSE') 
xlabel('Missing entries') 
hold on;
plot(Result_MC_mean(:,1),Result_MC_mean(:,2),'r-.','LineWidth',2);
if limY.active1
ylim([0,100])
end
hold off;
subplot(5,2,2)
plot(Result_MC(:,1),Result_MC(:,3),'b*');
title('Matrix Completion')
if limY.active1
ylim([0,100])
end
ylabel('NRMSE') 
xlabel('Missing entries') 
hold on;
plot(Result_MC_mean(:,1),Result_MC_mean(:,3),'r-.','LineWidth',2);
if limY.active1
ylim([0,100])
end
hold off;
end

if doLN
subplot(5,2,3)
plot(Result_Linear(:,1),Result_Linear(:,2),'b*');
title('Linear Interpolation')
if limY.active1
ylim([0,100])
end
ylabel('NMSE') 
xlabel('Missing entries') 
hold on;
plot(Result_Linear_mean(:,1),Result_Linear_mean(:,2),'r-.','LineWidth',2);
if limY.active1
ylim([0,100])
end
hold off;
subplot(5,2,4)
plot(Result_Linear(:,1),Result_Linear(:,3),'b*');
title('Linear Interpolation')
if limY.active1
ylim([0,100])
end
ylabel('NRMSE') 
xlabel('Missing entries') 
hold on;
plot(Result_Linear_mean(:,1),Result_Linear_mean(:,3),'r-.','LineWidth',2);
if limY.active1
ylim([0,100])
end
hold off;
end

if doSPL
subplot(5,2,5)
plot(Result_Spline(:,1),Result_Spline(:,2),'b*');
title('Spline Interpolation')
if limY.active1
ylim([0,100])
end
ylabel('NMSE') 
xlabel('Missing entries') 
hold on;
plot(Result_Spline_mean(:,1),Result_Spline_mean(:,2),'r-.','LineWidth',2);
if limY.active1
ylim([0,100])
end
hold off;
subplot(5,2,6)
plot(Result_Spline(:,1),Result_Spline(:,3),'b*');
title('Spline Interpolation')
if limY.active1
ylim([0,100])
end
ylabel('NRMSE') 
xlabel('Missing entries') 
hold on;
plot(Result_Spline_mean(:,1),Result_Spline_mean(:,3),'r-.','LineWidth',2);
if limY.active1
ylim([0,100])
end
hold off;
end

if doH
subplot(5,2,7)
plot(Result_Hermite(:,1),Result_Hermite(:,2),'b*');
title('Hermite Interpolation')
if limY.active1
ylim([0,100])
end
ylabel('NMSE') 
xlabel('Missing entries') 
hold on;
plot(Result_Hermite_mean(:,1),Result_Hermite_mean(:,2),'r-.','LineWidth',2);
if limY.active1
ylim([0,100])
end
hold off;
subplot(5,2,8)
plot(Result_Hermite(:,1),Result_Hermite(:,3),'b*');
title('Hermite Interpolation')
if limY.active1
ylim([0,100])
end
ylabel('NRMSE') 
xlabel('Missing entries') 
hold on;
plot(Result_Hermite_mean(:,1),Result_Hermite_mean(:,3),'r-.','LineWidth',2);
if limY.active1
ylim([0,100])
end
hold off;
end

if doMF
subplot(5,2,9)
plot(Result_missForest(:,1),Result_missForest(:,2),'b*');
title('missForest')
if limY.active1
ylim([0,100])
end
ylabel('NMSE') 
xlabel('Missing entries') 
hold on;
plot(Result_missForest_mean(:,1),Result_missForest_mean(:,2),'r-.','LineWidth',2);
if limY.active1
ylim([0,100])
end
hold off;
subplot(5,2,10)
plot(Result_missForest(:,1),Result_missForest(:,3),'b*');
title('missForest')
if limY.active1
ylim([0,100])
end
ylabel('NRMSE') 
xlabel('Missing entries') 
hold on;
plot(Result_missForest_mean(:,1),Result_missForest_mean(:,3),'r-.','LineWidth',2);
if limY.active1
ylim([0,100])
end
hold off;
end

print(['_Figures/dataset' idstr '_a_detailed.jpg'],'-djpeg','-r100')


%%













% figure(2)
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 16.0 8.0])
% lcell={};
% lcount=0;
% if doMC
% plot(Result_MC_mean(:,1),Result_MC_mean(:,2),line1.style,'LineWidth',line1.width);
% hold on;
% lcount=lcount+1;
% lcell{lcount}='MC';
% end
% if doGMC
% plot(Result_MC_L_mean(:,1),Result_MC_L_mean(:,2),line2.style,'LineWidth',line2.width);
% if limY.active2
% ylim([0,limY.value2])
% end
% hold on;
% lcount=lcount+1;
% lcell{lcount}='Laplacian Matrix Completion';
% end
% if doLN
% plot(Result_Linear_mean(:,1),Result_Linear_mean(:,2),'-.','LineWidth',2);
% if limY.active2
% ylim([0,limY.value2])
% end
% hold on;
% lcount=lcount+1;
% lcell{lcount}='Lin';
% end
% if doKNN
% plot(Result_KNN_mean(:,1),Result_KNN_mean(:,2),line3.style,'LineWidth',1);
% if limY.active2
% ylim([0,limY.value2])
% end
% hold on;
% lcount=lcount+1;
% lcell{lcount}='KNN';
% end
% if doSPL
% plot(Result_Spline_mean(:,1),Result_Spline_mean(:,2),'-*','LineWidth',1);
% if limY.active2
% ylim([0,limY.value2])
% end
% hold on;
% lcount=lcount+1;
% lcell{lcount}='Spline';
% end
% if doH
% plot(Result_Hermite_mean(:,1),Result_Hermite_mean(:,2),'-*','LineWidth',1);
% if limY.active2
% ylim([0,limY.value2])
% end
% hold on;
% lcount=lcount+1;
% lcell{lcount}='Hermite';
% end
% if doMF
% plot(Result_missForest_mean(:,1),Result_missForest_mean(:,2),line4.style,'LineWidth',1);
% lcount=lcount+1;
% lcell{lcount}='MF';
% end
% if limY.active2
% ylim([0,limY.value2])
% end
% %title('NMSE')
% ylabel('NMSE') 
% xlabel('Missing entries %') 
% grid on
% legend(lcell,'Location','northwest','Orientation','vertical')
% hold off;

if showPE
figure(4)
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4.0 4.0])
%subplot(3,1,2)
lcell={};
lcount=0;
if doMC
plot(Result_MC_mean(:,1),Result_MC_mean(:,4),'--','LineWidth',2);
hold on;
lcount=lcount+1;
lcell{lcount}='Matrix Completion';
end

if doGMC
plot(Result_MC_L_mean(:,1),Result_MC_L_mean(:,4),'-.','LineWidth',2);
hold on;
lcount=lcount+1;
lcell{lcount}='Laplacian Matrix Completion';
end

if doLN
plot(Result_Linear_mean(:,1),Result_Linear_mean(:,4),'-.','LineWidth',2);
hold on;
lcount=lcount+1;
lcell{lcount}='Lin';
end

if doKNN
plot(Result_KNN_mean(:,1),Result_KNN_mean(:,4),'-.','LineWidth',1);
if limY
ylim([0,0.01])
end
hold on;
lcount=lcount+1;
lcell{lcount}='KNN';
end

if doSPL
plot(Result_Spline_mean(:,1),Result_Spline_mean(:,4),'-*','LineWidth',1);
hold on;
lcount=lcount+1;
lcell{lcount}='Spline';
end

if doH
plot(Result_Hermite_mean(:,1),Result_Hermite_mean(:,4),'-*','LineWidth',1);
hold on;
lcount=lcount+1;
lcell{lcount}='Hermite';
end

if doMF
plot(Result_missForest_mean(:,1),Result_missForest_mean(:,4),'-*','LineWidth',1);
lcount=lcount+1;
lcell{lcount}='MF';
end
if limY
ylim([0,0.01])
end
title('PE')
ylabel('PE %') 
xlabel('Missing entries %') 
grid on
legend(lcell,'Location','northwest','Orientation','vertical')
hold off;
end

%%

%boxplot()

figure(2)
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 16.0 8.0])
%subplot(3,1,3)
lcell={};
lcount=0;
if doMC
plot(Result_MC_mean(:,1),Result_MC_mean(:,2),line1.style,'LineWidth',line1.width,'MarkerSize',14);
hold on;
lcount=lcount+1;
lcell{lcount}='Matrix Completion';
end
if doGMC
plot(Result_MC_L_mean(:,1),Result_MC_L_mean(:,2),line2.style,'LineWidth',line2.width,'MarkerSize',14);
if limY.active2
ylim([0,limY.value2])
end
hold on;
lcount=lcount+1;
lcell{lcount}='Laplacian Matrix Completion';
end
if doLN
plot(Result_Linear_mean(:,1),Result_Linear_mean(:,2),'-.','LineWidth',2);
if limY.active2
ylim([0,limY.value2])
end
hold on;
lcount=lcount+1;
lcell{lcount}='Lin';
end
if doKNN
plot(Result_KNN_mean(:,1),Result_KNN_mean(:,2),line3.style,'LineWidth',line3.width,'MarkerSize',14);
if limY.active2
ylim([0,limY.value2])
end
hold on;
lcount=lcount+1;
lcell{lcount}='KNN';
end
if doSPL
plot(Result_Spline_mean(:,1),Result_Spline_mean(:,2),'-*','LineWidth',1);
if limY.active2
ylim([0,limY.value2])
end
hold on;
lcount=lcount+1;
lcell{lcount}='Spline';
end
if doH
plot(Result_Hermite_mean(:,1),Result_Hermite_mean(:,2),'-*','LineWidth',1);
if limY.active2
ylim([0,limY.value2])
end
hold on;
lcount=lcount+1;
lcell{lcount}='Hermite';
end
if doMF
plot(Result_missForest_mean(:,1),Result_missForest_mean(:,2),line4.style,'LineWidth',line4.width,'MarkerSize',14);
lcount=lcount+1;
lcell{lcount}='MF';
end
if limY.active2
ylim([0,limY.value2])
end
%title('NMSE')
ylabel('NMSE') 
xlabel('Missing entries %') 
set(gca,'FontSize',36)
grid on
legend(lcell,'Location','northwest','Orientation','vertical','FontSize',36)
hold off;
print(['_Figures/' idstr '_NMSE.png'],'-dpng','-r100')

%%

figure(3)
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 16.0 9.0])
%subplot(3,1,2)
lcell={};
lcount=0;
if doMC
plot(Result_MC_mean(:,1),Result_MC_mean(:,3),line1.style,'LineWidth',line1.width,'MarkerSize',14);
hold on;
lcount=lcount+1;
lcell{lcount}='Matrix Completion';
end

if doGMC
plot(Result_MC_L_mean(:,1),Result_MC_L_mean(:,3),line2.style,'LineWidth',line2.width,'MarkerSize',14);
hold on;
lcount=lcount+1;
lcell{lcount}='Laplacian Matrix Completion';
end

if doLN
plot(Result_Linear_mean(:,1),Result_Linear_mean(:,3),'-.','LineWidth',2);
hold on;
lcount=lcount+1;
lcell{lcount}='Lin';
end

if doKNN
plot(Result_KNN_mean(:,1),Result_KNN_mean(:,3),line3.style,'LineWidth',line3.width,'MarkerSize',14);
if limY.active2
ylim([0,limY.value2])
end
hold on;
lcount=lcount+1;
lcell{lcount}='KNN';
end

if doSPL
plot(Result_Spline_mean(:,1),Result_Spline_mean(:,3),'-*','LineWidth',1);
hold on;
lcount=lcount+1;
lcell{lcount}='Spline';
end

if doH
plot(Result_Hermite_mean(:,1),Result_Hermite_mean(:,3),'-*','LineWidth',1);
hold on;
lcount=lcount+1;
lcell{lcount}='Hermite';
end

if doMF
plot(Result_missForest_mean(:,1),Result_missForest_mean(:,3),line4.style,'LineWidth',line4.width,'MarkerSize',14);
lcount=lcount+1;
lcell{lcount}='MF';
end
if limY.active2
ylim([0,limY.value2])
end
%title('NRMSE')
ylabel('NRMSE') 
xlabel('Missing entries %') 
set(gca,'FontSize',36)
grid on
legend(lcell,'Location','northwest','Orientation','vertical','FontSize',36)
hold off;
print(['_Figures/' idstr '_NRMSE.png'],'-dpng','-r100')
%%

if showHist
figure(5)
%set(gcf,'PaperUnits','inches','PaperPosition',[0 0 4.0 4.0])
%subplot(3,1,1)
lcell={};
lcount=0;
if doMC
plot(Result_MC_mean(:,1),Result_MC_mean(:,5),'--','LineWidth',2);
hold on;
lcount=lcount+1;
lcell{lcount}='Matrix Completion';
end

if doGMC
plot(Result_MC_L_mean(:,1),Result_MC_L_mean(:,5),'-.','LineWidth',2);
if limY.active2
ylim([0,limY.value2])
end
hold on;
lcount=lcount+1;
lcell{lcount}='Laplacian matrix completion';
end

if doLN
plot(Result_Linear_mean(:,1),Result_Linear_mean(:,5),'-.','LineWidth',2);
if limY.active2
ylim([0,limY.value2])
end
hold on;
lcount=lcount+1;
lcell{lcount}='Lin';
end

if doKNN
plot(Result_KNN_mean(:,1),Result_KNN_mean(:,5),'-.','LineWidth',1);
if limY.active2
ylim([0,limY.value2])
end
hold on;
lcount=lcount+1;
lcell{lcount}='KNN';
end

if doSPL
plot(Result_Spline_mean(:,1),Result_Spline_mean(:,5),'-*','LineWidth',1);
if limY.active2
ylim([0,limY.value2])
end
hold on;
lcount=lcount+1;
lcell{lcount}='Spline';
end

if doH
plot(Result_Hermite_mean(:,1),Result_Hermite_mean(:,5),'-*','LineWidth',1);
if limY.active2
ylim([0,limY.value2])
end
hold on;
lcount=lcount+1;
lcell{lcount}='Hermite';
end

if doMF
plot(Result_missForest_mean(:,1),Result_missForest_mean(:,5),'-*','LineWidth',1);
lcount=lcount+1;
lcell{lcount}='MF';
end

if limY.active2
ylim([0,limY.value2])
end
title('Histograms Percentage Error')
ylabel('PE %') 
xlabel('Missing entries %') 
grid on
legend(lcell,'Location','northwest','Orientation','vertical','FontSize',12)
hold off;

print(['_Figures/dataset' idstr '_b_comparative.jpg'],'-djpeg','-r100')
end

%% Some conclusions

% M=Result_Linear(:,2)-Result_MC(:,2)
% Mind=find(M>0)


% 
col=3;
figure(8)
if doMC
boxplot([...
Result_MC(Result_MC(:,1)==5,col)...
,Result_MC(Result_MC(:,1)==10,col)...
,Result_MC(Result_MC(:,1)==15,col)...
,Result_MC(Result_MC(:,1)==20,col)...
,Result_MC(Result_MC(:,1)==25,col)...
,Result_MC(Result_MC(:,1)==30,col)...
,Result_MC(Result_MC(:,1)==35,col)...
,Result_MC(Result_MC(:,1)==40,col)...
% ,Result_MC(Result_MC(:,1)==45,col)...
% ,Result_MC(Result_MC(:,1)==50,col)...
],'Colors' ,'g','BoxStyle','filled','Labels',{
'5'...
,'10'...
,'15'...
,'20'...
,'25'...
,'30'...
,'35'...
,'40'...
% ,'45'...
% ,'50'...
})
xlabel('Missing entries %','FontSize', 12) 
ylabel('NRMSE','FontSize', 12) 
 box_vars = findall(gca,'Tag','Box');
hold on;
end

if doMF
boxplot([...
    Result_missForest(Result_missForest(:,1)==5,col)...
    ,Result_missForest(Result_missForest(:,1)==10,col)...
    ,Result_missForest(Result_missForest(:,1)==15,col)...
    ,Result_missForest(Result_missForest(:,1)==20,col)...
    ,Result_missForest(Result_missForest(:,1)==25,col)...
    ,Result_missForest(Result_missForest(:,1)==30,col)...
    ,Result_missForest(Result_missForest(:,1)==35,col)...
    ,Result_missForest(Result_missForest(:,1)==40,col)...
%     ,Result_missForest(Result_missForest(:,1)==45,col)...
%     ,Result_missForest(Result_missForest(:,1)==50,col)...
],'Colors' ,'y','BoxStyle','filled','Labels',{
'5'...
,'10'...
,'15'...
,'20'...
,'25'...
,'30'...
,'35'...
,'40'...
% ,'45'...
% ,'50'...
})
xlabel('Missing entries %','FontSize', 12) 
ylabel('NRMSE','FontSize', 12) 
 box_vars = findall(gca,'Tag','Box');
hold on;
end

if doGMC
boxplot([...
Result_MC_L(Result_MC_L(:,1)==5,col)...
,Result_MC_L(Result_MC_L(:,1)==10,col)...
,Result_MC_L(Result_MC_L(:,1)==15,col)...
,Result_MC_L(Result_MC_L(:,1)==20,col)...
,Result_MC_L(Result_MC_L(:,1)==25,col)...
,Result_MC_L(Result_MC_L(:,1)==30,col)...
,Result_MC_L(Result_MC_L(:,1)==35,col)...
,Result_MC_L(Result_MC_L(:,1)==40,col)...
% ,Result_MC_L(Result_MC_L(:,1)==45,col)...
% ,Result_MC_L(Result_MC_L(:,1)==50,col)...
],'Colors' ,'m','BoxStyle','filled','Labels',{
'5'...
,'10'...
,'15'...
,'20'...
,'25'...
,'30'...
,'35'...
,'40'...
% ,'45'...
% ,'50'...
})
xlabel('Missing entries %','FontSize', 12) 
ylabel('NRMSE','FontSize', 12) 
 box_vars = findall(gca,'Tag','Box');
hold on;
end


if doKNN
boxplot([Result_KNN(Result_KNN(:,1)==5,col)...
    ,Result_KNN(Result_KNN(:,1)==10,col)...
    ,Result_KNN(Result_KNN(:,1)==15,col)...
    ,Result_KNN(Result_KNN(:,1)==20,col)...
    ,Result_KNN(Result_KNN(:,1)==25,col)...
    ,Result_KNN(Result_KNN(:,1)==30,col)...
    ,Result_KNN(Result_KNN(:,1)==35,col)...
    ,Result_KNN(Result_KNN(:,1)==40,col)...
%     ,Result_KNN(Result_KNN(:,1)==45,col)...
%     ,Result_KNN(Result_KNN(:,1)==50,col)...
],'Colors' ,'c','BoxStyle','filled','Labels',{
'5'...
,'10'...
,'15'...
,'20'...
,'25'...
,'30'...
,'35'...
,'40'...
% ,'45'...
% ,'50'...
})
xlabel('Missing entries %','FontSize', 12) 
ylabel('NRMSE','FontSize', 12) 
 box_vars = findall(gca,'Tag','Box');
legend(box_vars([1,11,21,31]), {'KNN','Laplacian matrix completion','MissForest','Matrix completion'},'Location','NorthWest','FontSize',12);
hold on;
end

print(['_Figures/' idstr '_NRMSE_boxplot.png'],'-dpng','-r100')
