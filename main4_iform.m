% Generate environmental contour lines using iForm with IEC and gmm methods
clear;
clc;
close all;
addpath('../gmm');
%% load data x1 and x2
% data_filename = 'onshore_detrend_45_135.mat';  % Data file to load
data_filename = 'offshore_detrend.mat';  % Data file to load
n_points = 100; 

f = fullfile('../data/hovsore_wind/data',data_filename);
calmstd = load(f).calmstd;
x1 = calmstd.u;
x2 = calmstd.u_stdv;
% calculate z1, z2, which are normal space correlated variables from iForm
n_m = 50*365*24*6;
pd = makedist('Normal','mu',0,'sigma',1);
betav = icdf(pd, 1-1/n_m);
z1 = linspace(-3,4.5, n_points);
z2_p = sqrt(betav.^2-z1.^2);
z2_n = -sqrt(betav.^2-z1.^2);

%% IEC with Iref = 0.12, Weibull and Lognormal distribution
Iref = 0.12;
pd_f1 = fitdist(x1, 'Weibull'); %fitted pd1
u = icdf('Weibull',normcdf(z1),pd_f1.A,pd_f1.B);
m = Iref.*(0.75.*u+3.8);
v = (1.4*Iref).^2;
mu_ln = log((m.^2)./sqrt(v+m.^2));
sigma_ln = sqrt(log(v./(m.^2)+1));
sigma1 = icdf('Lognormal',normcdf(z2_p),mu_ln,sigma_ln);
%negative
sigma2 = icdf('Lognormal',normcdf(z2_n),mu_ln,sigma_ln);
filename = sprintf('res/contour_iec_iref_%.2f_%s', Iref, data_filename);
save(filename, "u", "sigma2", "sigma1");


%% IEC with Iref from wind measurement, Weibull and Lognormal distribution
ti = x2./x1;
index = (x1>=14) & (x1<16);
sig_t = x2(index);
ti_t = ti(index);
Iref = mean(ti_t);
m = Iref.*(0.75.*u+3.8);
v = (1.4*Iref).^2;
mu_ln = log((m.^2)./sqrt(v+m.^2));
sigma_ln = sqrt(log(v./(m.^2)+1));
sigma1 = icdf('Lognormal',normcdf(z2_p),mu_ln,sigma_ln);
%negative
sigma2 = icdf('Lognormal',normcdf(z2_n),mu_ln,sigma_ln);
filename = sprintf('res/contour_iec_iref_%.2f_%s', Iref, data_filename);
save(filename, "u", "sigma2", "sigma1");



%% GMM
f = sprintf('res/gmm_model_%s', data_filename);
gmModel = load(f).gmModel;
pd = gmm_df(gmModel);
mxpdf = pd.mxpdf;
mxpdf2 = pd.mxpdf2;
pdfx1 = pd.pdfx1;
pdfx2 = pd.pdfx2;
cdfx1 = pd.cdfx1;
cdfx2 = pd.cdfx2;

% get inverse cdf of x1, thus get mean wind speed u
n = 1e5;
x1t = linspace(0, max(x1)+3*std(x1), n);
% x2t = linspace(0, max(x2)+3*std(x2), n);
y1t = cdfx1(x1t);
% y2t = cdfx2(x2t);
icdfx1 = @(x) interp1(y1t,x1t,x);
% icdfx2 = @(x) interp1(y2t,x2t,x);
u = icdfx1(normcdf(z1));

sigma1 =zeros(1, n_points);
for k = 1:numel(u)
    ut = u(k);
    pdfc = @(x) mxpdf2(ut,x)./pdfx1(ut);
    cdfc = @(x) integral(pdfc,-5,x);
    sigma_lb = 0;
    sigma_ub = 5;
    t = normcdf(z2_p(k)); 
    % Define the function handle for g(x) = cdfc(x) - t
    % find the root x, so that cdfc(x) = t;
    g = @(x) cdfc(x) - t;
    % Call the general bisection method
    tolerance = 1e-6;
    sigma1(k) = bisection_method(sigma_lb, sigma_ub, g, tolerance);  
end

%negative

sigma2 =zeros(1, numel(z1));
for k = 1:numel(u)
    ut = u(k);
    pdfc = @(x) mxpdf2(ut,x)./pdfx1(ut);
    cdfc = @(x) integral(pdfc,-3,x);
    sigma_lb = 0;
    sigma_ub = 3;
    t = normcdf(z2_n(k));
    % Define the function handle for g(x) = fun(x) - t
    g = @(x) cdfc(x) - t;
    % Call the general bisection method
    tolerance = 1e-6;
    sigma2(k) = bisection_method(sigma_lb, sigma_ub, g, tolerance);        
end

filename = sprintf('res/contour_gmm_%s', data_filename);
save(filename, "u", "sigma2", "sigma1");







