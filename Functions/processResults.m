
close all;
R=csvread('CSV/nmse_vs_missing_uniform.csv');

plot(R(:,1),R(:,2),'-o','LineWidth',2,'MarkerSize',4);
ylim([0 100]);
xlabel({'Missing entries uniformly distributed(%) '})
ylabel({'NMSE(%)'})
title('NMSE vs Missing entries for a single day')