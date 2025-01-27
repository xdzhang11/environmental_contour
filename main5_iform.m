%% IFORM analysis
% Reliability index beta
clear;
clc;
calmstd = load('../data/hovsore_wind/data/offshore_detrend.mat').calmstd;
x1 = calmstd.u;
x2 = calmstd.u_stdv;

figure;
plot(x1,x2,'.');
xlabel('X_1')
ylabel('X_2')
ylim([0 7])
xlim([0 40])
title('50-year turbulence')

hold on 

n_m = 50*365*24*6;
mu = 0;
sigma = 1;
pd = makedist('Normal','mu',mu,'sigma',sigma);
beta = icdf(pd, 1-1/n_m);

u1 = normrnd(0,1, [1e5,1]);
u2 = sqrt(beta.^2-u1.^2);

pd_f1 = fitdist(x1, 'Weibull'); %fitted pd1
u = icdf('Weibull',normcdf(u1),pd_f1.A,pd_f1.B);

% ti = x2./x1;
% index = (x1>=14) & (x1<16);
% sig_t = x2(index);
% ti_t = ti(index);
% Iref = mean(ti_t);
Iref = 0.12;

m = Iref.*(0.75.*u+3.8);
std = 1.4*Iref;
v = std.^2;
mu_ln = log((m.^2)./sqrt(v+m.^2));
sigma_ln = sqrt(log(v./(m.^2)+1));
sigma1 = icdf('Lognormal',normcdf(u2),mu_ln,sigma_ln);

scatter(u,sigma1,'b+')
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
%negative
u2 = -sqrt(beta.^2-u1.^2);

sigma2 = icdf('Lognormal',normcdf(u2),mu_ln,sigma_ln);
hold on
scatter(u,sigma2,'b+')
xlabel('U (m/s)')
ylabel('\sigma_u (m/s)')

u2 = sqrt(beta.^2-u1.^2);
ti = x2./x1;
index = (x1>=14) & (x1<16);
sig_t = x2(index);
ti_t = ti(index);
Iref = mean(ti_t);

m = Iref.*(0.75.*u+3.8);
std = 1.4*Iref;
v = std.^2;
mu_ln = log((m.^2)./sqrt(v+m.^2));
sigma_ln = sqrt(log(v./(m.^2)+1));
sigma1 = icdf('Lognormal',normcdf(u2),mu_ln,sigma_ln);

scatter(u,sigma1,'y+')
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
%negative
u2 = -sqrt(beta.^2-u1.^2);

sigma2 = icdf('Lognormal',normcdf(u2),mu_ln,sigma_ln);
hold on
scatter(u,sigma2,'y+')
xlabel('U (m/s)')
ylabel('\sigma_u (m/s)')



