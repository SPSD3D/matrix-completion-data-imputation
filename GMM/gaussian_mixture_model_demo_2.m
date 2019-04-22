% Test out GAUSSIAN_MIXTURE_MODEL
%
% file:      	gaussian_mixture_model_test.m, (c) Matthew Roughan, Tue Jul 21 2009
% directory:   /home/mroughan/src/matlab/NUMERICAL_ROUTINES/
% created: 	Tue Jul 21 2009 
% author:  	Matthew Roughan 
% email:   	matthew.roughan@adelaide.edu.au
% 
%
%
clear;

% set parameters for 2 classes, similar to Internet traffic data
C = 2;
w = [0.7 0.3];
mu = [1 4];
sigma = [0.9 1.9];
x = -4:0.25:14;
cdf = w(1) * normcdf(x, mu(1), sigma(1)) + ...
      w(2) * normcdf(x, mu(2), sigma(2));

% generate a set of data with 2 classes
N = 10000;
r = 1 + (rand(N,1) > w(1)); % which mixture should it come from
mu_d = mu(r);
sigma_d = sigma(r);
values = mu_d + sigma_d .* randn(1,N);
plot(values)




id=2;
idstr = sprintf('%03d',id);
C = 3;
datasetfolder='CSV/';
dataset=['ProcessedDataset' idstr '/'];
resultsfolder='CSV/Results/';
datasettarget=['dataset' idstr];

writecompleted=true;
completedfolder=[datasetfolder dataset 'Completed/'];

Data=csvread([datasetfolder dataset 'initial.csv']);



x=Data(:);
Fs=2/60
N = length(x);
xdft = fft(x);
xdft = xdft(1:N/2+1);
psdx = (1/(2*pi*N)) * abs(xdft).^2;
psdx =  abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:(2*pi)/N:pi;
psd=10*log10(psdx)
psd=psd/sum(psd)
freq=freq/pi
plot(freq,psd,'o')
grid on
title('Periodogram Using FFT')
xlabel('Normalized Frequency (\times\pi rad/sample)') 
ylabel('Power/Frequency (dB/rad/sample)')




x=0:.000001:0.01;
values=psd'



C=3




% now estimate the parametes of the distribution
[mu_est, sigma_est, w_est, counter, difference] = gaussian_mixture_model(values, C, 1.0e-3);
mu_est'
sigma_est'
w_est'

% compare empirical data to estimated distribution
p1_est = w_est(1) * norm_density(x, mu_est(1), sigma_est(1));
p2_est = w_est(2) * norm_density(x, mu_est(2), sigma_est(2));
p3_est = w_est(3) * norm_density(x, mu_est(3), sigma_est(3));


figure(3)
plot2 = plot(x, p1_est+p2_est+p3_est, 'r--', 'linewidth', 2);
set(gca,'linewidth', 2, 'fontsize', 16)



