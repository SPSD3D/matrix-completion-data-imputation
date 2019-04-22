function [fin,p1_est,p2_est,p3_est,x] = histgmmv2( Din )
%HISTGMM Summary of this function goes here
%   Detailed explanation goes here
values=Din(:);
values=values';
x=min(values):max(values);
C=3;

GMModel = fitgmdist(values',3)
mu_est=GMModel.mu;
w_est=GMModel.ComponentProportion;
sigma_est=GMModel.Sigma(:)
% compare empirical data to estimated distribution
p1_est = w_est(1) * norm_density(x, mu_est(1), sigma_est(1));
p2_est = w_est(2) * norm_density(x, mu_est(2), sigma_est(2));
p3_est = w_est(3) * norm_density(x, mu_est(3), sigma_est(3));
% p1 = w(1) * norm_density(x, mu(1), sigma(1));
% p2 = w(2) * norm_density(x, mu(2), sigma(2));

fin=p1_est+p2_est+p3_est

end

