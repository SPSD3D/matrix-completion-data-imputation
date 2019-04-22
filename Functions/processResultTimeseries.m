MI=csvread('CSV/MInitial.csv');

MC=csvread('CSV/MCompletion.csv');

MI2=reshape(MI,1,[]);
MC=reshape(MC,1,[]);

figure(4)
plot(MI2,'--','LineWidth',3)
hold on;
plot(MC,'-o','LineWidth',1,'MarkerSize',2)
hold off;
legend({'Initial timeseries','Matrix completion'},'FontSize',14)

