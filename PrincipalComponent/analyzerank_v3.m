clc;
clear variables;
close all;

datatimeseriesfolder=['CSV/'];
dataset='ProcessedDataset001/';
file='initial.csv';
Din=csvread([datatimeseriesfolder dataset file]);



[fin,p1_est,p2_est,p3_est,x]=histgmm(Din);






Rhist=histogram(Din(:));
figure(3)
Rhistval=Rhist.Values/sum(Rhist.Values)
plot(x,Rhistval)
hold on;


plot1 = plot(x, p1_est, 'b-', 'linewidth', 2);
hold on
plot2 = plot(x, p2_est, 'g-.', 'linewidth', 2);
hold on
plot3 = plot(x, p3_est, 'm-.', 'linewidth', 2);
hold on
plot4 = plot(x, fin, 'r--', 'linewidth', 2);



%set(gcf, 'PaperPosition', [0 0 5 3]);
%print('-dpng', 'gaussian_mixture_model.png');


