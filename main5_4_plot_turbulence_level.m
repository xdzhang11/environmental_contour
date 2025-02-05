% Plot turbulence level at different wind speeds
clear;
clc;
close all; 

% load data x1 and x2
filename = 'offshore_detrend.mat';
f = fullfile('../../../data/hovsore_wind/data',filename);
calmstd = load(f).calmstd;
u = calmstd.u;
stdv = calmstd.u_stdv;

numInterval = 11;

sig_ap = zeros(1,numInterval-1);
sig_a = zeros(1,numInterval-1);
sig_b = zeros(1,numInterval-1);
sig_c = zeros(1,numInterval-1);
sig_IEC = zeros(1,numInterval-1);
% 90 quantile of standard deviation at wind speed 15m/s
ti = stdv./u;
index = (u>=14) & (u<16);
sig_t = stdv(index);
ti_t = ti(index);
Iref = quantile(ti_t, 0.7);

for k = 1:numInterval-1
    u1 = 2*k+1;
    u2 = u1+2;
    index = (u>=u1) & (u<u2);
    sig_t = stdv(index);
    sig_90q(k) = quantile(sig_t, 0.9);
    v = 0.5*(u1+u2); %wind speed
    sig_IEC(k) = Iref*(0.75*v+5.6);
    sig_ap(k) = 0.18*(0.75*v+5.6);
    sig_a(k) = 0.16*(0.75*v+5.6);
    sig_b(k) = 0.14*(0.75*v+5.6);
    sig_c(k) = 0.12*(0.75*v+5.6);
    x_p(k)=v;
end
figure;
hold on
scatter(x_p, sig_90q, 'DisplayName', 'Representative (90% quantile)');
plot(x_p, sig_IEC, '-.' ,'DisplayName', 'IEC (Iref = 0.066, 70% quantile at 15m/s)');
plot(x_p, sig_ap, 'DisplayName', 'A+ (Iref = 0.18)');
plot(x_p, sig_a, 'DisplayName', 'A (Iref = 0.16)');
plot(x_p, sig_b, 'DisplayName', 'B (Iref = 0.14)');
plot(x_p, sig_c, 'DisplayName', 'C (Iref = 0.12)');
xlabel('\it{u} (m/s)')
ylabel('\sigma_1 (m/s)')
box on;
legend( 'Location','northwest' );