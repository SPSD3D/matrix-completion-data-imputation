clc;
clear;
close all;

datatimeseriesfolder=['CSV/'];
dataset='ProcessedDataset002/';
file='initial.csv';
randomfile='uniform_missing_30_0002.csv';
%randomfile='perblock_missing_30_0002.csv';
D=csvread([datatimeseriesfolder dataset file]);
numel(D)


[L,Sm]=RobustPCA(D)
pErr=100*norm(D-L)/norm(D)
dHist=100*(mean((D(:)-L(:)).^2)/mean((D(:)-mean(D(:))).^2));




% 
% 
% 
% D=D(:)
% D=D(1:64)
% 
% b=hnk(D,32);
% [L,Sm]=RobustPCA(b)
% 
% figure(4)
% [U,S,V]=svd(b);
% SD=diag(S);
% plot(SD/sum(SD),'LineWidth',3);
% grid on
% 
% figure(5)
% [U,S,V]=svd(L-mean(mean(L)));
% SD=diag(S);
% plot(SD/sum(SD),'LineWidth',3);
% grid on
% 
% E=csvread([datatimeseriesfolder dataset randomfile]);
% figure(6)
% % imagesc(E);
% H=im2bw(E, 0.01);
% imagesc(H)
