clc;
clear;
close all;

M8=csvread('MInitialTotalMatrixTrunc8ACD.csv');
[U,S,V]=svd(M8);
SD=diag(S);
SF8=SD/sum(SD);
figure(8)
plot(SF8,'LineWidth',2)
set(gca,'fontsize',14)


M7=csvread('MInitialTotalMatrixTrunc4FEV.csv');
[U,S,V]=svd(M7);
SD=diag(S);
SF7=SD/sum(SD);
figure(7)
plot(SF7,'LineWidth',2)
set(gca,'fontsize',14)




M3=csvread('MInitialSingleACD.csv');
[U,S,V]=svd(M3);
SD=diag(S);
SF3=SD/sum(SD);
figure(3)
plot(SF3,'LineWidth',2)
set(gca,'fontsize',14)

M4=csvread('MInitialSingleFEV.csv');
[U,S,V]=svd(M4);
SD=diag(S);
SF4=SD/sum(SD);
figure(4)
plot(SF4,'LineWidth',2)
set(gca,'fontsize',14)

M1=csvread('MInitialSimilarACD.csv');
[U,S,V]=svd(M1);
SD=diag(S);
SF1=SD/sum(SD);
figure(1)
plot(SF1,'LineWidth',2)
hold on
plot(SF3,'--','LineWidth',2)
hold off
legend({'Similar patients','Single patient'},'FontSize',14)
set(gca,'fontsize',14)


M2=csvread('MInitialSimilarFEV.csv');
[U,S,V]=svd(M2);
SD=diag(S);
SF2=SD/sum(SD);
figure(2)
plot(SF2,'LineWidth',2)




M5=csvread('MInitialTotalMatrixACD.csv');
[U,S,V]=svd(M5);
SD=diag(S);
SF5=SD/sum(SD);
figure(5)
plot(SF5,'LineWidth',2)
set(gca,'fontsize',14)

M6=csvread('MInitialTotalMatrixFEV.csv');
[U,S,V]=svd(M6);
SD=diag(S);
SF6=SD/sum(SD);
figure(6)
plot(SF6,'LineWidth',2)
set(gca,'fontsize',14)
