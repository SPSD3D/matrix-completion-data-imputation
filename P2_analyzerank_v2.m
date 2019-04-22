clc;
clear;
close all;

datasetid=50;

datatimeseriesfolder=['CSV/'];
dataset=['ProcessedDataset' sprintf('%03d',datasetid) '/' ];
file='initial.csv';
Din=csvread([datatimeseriesfolder dataset file]);
%Din=Din-mean(mean(Din));
D=Din;
numel(D);



% subplot(1,2,1)
[U,Sd,V]=svd(D);
SD=diag(Sd);
SDP=SD/sum(SD);
figure(3)
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 16 10])
bar(SDP);
xlabel('Component') % x-axis label
ylabel('Singular value') % y-axis label
set(gca,'FontSize',36)
%plot(SDP,'LineWidth',3);
grid on
print(['Figures/' num2str(datasetid) '-singular.png'],'-dpng','-r100')


SDPC=[];
for i=1:length(SDP)
    SDPC(i)=sum(SDP(1:i));
end
%subplot(1,2,2)
figure(4)
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 16 10])
bar(SDPC);
xlabel('Number of components') % x-axis label
ylabel('Sum of singular values') % y-axis label
set(gca,'FontSize',36)
ylim([0 1])
%plot(SDPC,'LineWidth',3);
grid on
set(findall(gcf,'-property','FontSize'),'FontSize',36)
print(['Figures/' num2str(datasetid) '-singular-sum.png'],'-dpng','-r100')

% figure(6)
% X=max(size(D,1),size(D,2));
% Xopt=X;
% Opt=D;
% Circ=0;
% for j=1:numel(D)
% D=circularshift(D);
% [U,Sd,V]=svd(D);
% SD=diag(Sd);
% SDP=SD/sum(SD);
% SDPC=[];
% for i=1:length(SDP)
%     SDPC(i)=sum(SDP(1:i));
% end
% Xcur=numel(find(SDPC(SDPC<0.85)));
% 
% if Xcur<Xopt
%     Xopt=Xcur;
%     Opt=D;
%     Circ=j;
% end
% end
% 
% [U,Sd,V]=svd(Opt);
% SD=diag(Sd);
% SDP=SD/sum(SD);
% plot(SDP,'LineWidth',3);
% grid on
% hold on
% [U,Sd,V]=svd(Din);
% SD=diag(Sd);
% SDP=SD/sum(SD);
% plot(SDP,'LineWidth',2);
% grid on
% hold off
