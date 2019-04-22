clc;
clear;
close all;

id=2;
idstr = sprintf('%03d',id);

dataset=['ProcessedDataset' idstr '/'];
datatimeseriesfolder=['CSV/' dataset];


figure(1)

j=150;

H=dir([datatimeseriesfolder '*.csv']);

init=csvread([datatimeseriesfolder H(1).name]);
init=reshape(init,1,[]);

missing=csvread([datatimeseriesfolder H(j).name]);
missing=reshape(missing,1,[]);
inds=find(missing);

mm=init;
mm(inds)=NaN;


% completed_linear=csvread([datatimeseriesfolder 'Completed/completed_linear_' H(j).name]);
% completed_linear=reshape(completed_linear,1,[]);
% 
% completed_spline=csvread([datatimeseriesfolder 'Completed/completed_spline_' H(j).name]);
% completed_spline=reshape(completed_spline,1,[]);
% 
% completed_hermite=csvread([datatimeseriesfolder 'Completed/completed_hermite_' H(j).name]);
% completed_hermite=reshape(completed_hermite,1,[]);

completed_mc=csvread([datatimeseriesfolder 'Completed/completed_mc_' H(j).name]);
completed_mc=reshape(completed_mc,1,[]);
% 
completed_mc_laplacian=csvread([datatimeseriesfolder 'Completed/completed_mc_laplacian' H(j).name]);
completed_mc_laplacian=reshape(completed_mc_laplacian,1,[]);
%
completed_missforest=csvread([datatimeseriesfolder 'Completed/completed_missforest_' H(j).name]);
completed_missforest=reshape(completed_missforest,1,[]);


completed_knn=csvread([datatimeseriesfolder 'Completed/completed_KNN_' H(j).name]);
completed_knn=reshape(completed_knn,1,[]);

plot(mm,'mo',...
    'LineWidth',2,...
    'MarkerEdgeColor','k',...
    'MarkerFaceColor',[.49 1 .63],...
    'MarkerSize',10)
hold on;

% plot(completed_linear,'LineWidth',2)
% hold on;
% plot(completed_spline,'LineWidth',2)
% hold on;
% plot(completed_hermite,'LineWidth',2)
% hold on;

plot(completed_mc,'-.','LineWidth',2)
hold on;
plot(completed_mc_laplacian,'-.','LineWidth',2)
hold on;
plot(completed_missforest,'-.','LineWidth',2)
hold on;
plot(completed_knn,'-.','LineWidth',2)
hold on;
plot(init,'-','LineWidth',2)
hold on;
% plot(completed_missforest)
legend({'Missing entries','Matrix Completion','Laplacian Matrix Completion','MissForest','KNN','Ground Truth','Linear interpolation','Spline fit', 'Hermite interpolation','Missforest'},'Location','NorthEast','FontSize',18);
hold off;

