% Generate correlated sample from Nataf and copulas
clear;
clc;
addpath('nataf');
addpath('../gmm');
addpath('../utilities');

% load meaurement data x1, x2

data_filename = 'offshore_detrend.mat';  % Data file to load
f = fullfile('../data/hovsore_wind/data',data_filename);
calmstd = load(f).calmstd;
x1 = calmstd.u;
x2 = calmstd.u_stdv;


n_sample = numel(x1); % sample size

mu_i = mean(x1);
sigma_i = std(x1);
mu_j = mean(x2);
sigma_j = std(x2);

rho_x_m = corrcoef(x1, x2);
rho_x = rho_x_m(1,2);
fprintf('The correlation of coefficient of x1 and x2 is %2f\n', rho_x);

d1.name = 'Weibull';
pd_f1 = fitdist(x1, 'Weibull'); %fitted pd1
d1.parameter = [pd_f1.A,pd_f1.B];
d2.name = 'Lognormal';
ln_mu_j = log(mu_j^2/sqrt(sigma_j^2+mu_j^2));
ln_sigma_j = sqrt(log(sigma_j^2/mu_j^2+1));
d2.parameter = [ln_mu_j, ln_sigma_j];
x_stat = [mu_i, mu_j, sigma_i, sigma_j];

%% Nataf
d1 = @(x) icdf(d1.name, x, d1.parameter(1), d1.parameter(2)); % inverse cdf function for x1
d2 = @(x) icdf(d2.name, x, d2.parameter(1), d2.parameter(2)); % inverse cdf function for x2

rho_z = nataffit(rho_x, d1, d2, x_stat);

xt = nattafsample(rho_z, d1, d2, n_sample);
x_i = xt(:,1);
x_j = xt(:,2);

rho_x_MCS = corrcoef(x_i, x_j);
rho_x_MCS = rho_x_MCS(1,2);
fprintf('The correlation of coefficient estimated from nataf MCS sample is %2f\n', rho_x_MCS);

filename = sprintf('res/sample_nataf_%s.mat', data_filename);
save(filename, 'x_i','x_j')

%% Copula
u = cdf('Weibull',x1,pd_f1.A,pd_f1.B);
v = cdf('Lognormal',x2,ln_mu_j,ln_sigma_j);

% Gaussian copula
% rhohat = copulafit('Gaussian',[u v]);
% r = copularnd('Gaussian',rhohat,n_sample);
% file_name = 'res/sample_copula_gaussian.mat';
%--------------------------------------------------------------------------
% t copula
% [Rho,nu] = copulafit('t',[u v]);
% r = copularnd('t',Rho,nu,n_sample);
% file_name = 'res/sample_copula_t.mat';
%--------------------------------------------------------------------------
% Gumbel copula
paramhat = copulafit('Gumbel',[u v]);
r = copularnd('Gumbel',paramhat,n_sample); % correlated sample
filename = sprintf('res/sample_copula_gumbel_%s.mat', data_filename);

u_c = r(:,1);
v_c = r(:,2);

x_i = icdf('Weibull',u_c,pd_f1.A,pd_f1.B);
x_j = icdf('Lognormal',v_c,ln_mu_j,ln_sigma_j);
rho_x_m = corrcoef(x_i, x_j);
rho_x_MCS = rho_x_m(1,2);

fprintf('The correlation of coefficient estimated from copula MCS sample is %2f\n', rho_x_MCS);
save(filename, 'x_i','x_j')

%% IEC sample
x_i = wblrnd(pd_f1.A,pd_f1.B, [n_sample,1]);

Iref = 0.0745;

a = Iref.*(0.75.*x_i+3.3);
b = 0.27.*x_i+1.4;

x_j = wblrnd(a , b);

roh_x_MCS = (mean(x_i.*x_j)-mu_i*mu_j)/sigma_i/sigma_j;
fprintf('The correlation of coefficient estimated from IEC MCS sample is %2f\n', roh_x_MCS);
filename = sprintf('res/sample_iec_%s.mat', data_filename);
save(filename, 'x_i','x_j')


%% GMM sample

model_filename = sprintf('gmm_model_%s', data_filename);  % Modify filename for saving model

f = fullfile('res', model_filename);
gmModel = load(f).gmModel;
X_p = gmm_sample(gmModel, n_sample);
x_i = X_p(:,1);
x_j = X_p(:,2);
filename = sprintf('res/sample_gmm_%s.mat', data_filename);
save(filename, 'x_i','x_j')
