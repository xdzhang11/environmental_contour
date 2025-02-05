% Plot environmental contours
clear;
clc;
close all;
% figure;
% data_filename = 'onshore_detrend_45_135.mat';  % Data file to load
data_filename = 'offshore_detrend.mat';  % Data file to load
f = fullfile('../data/hovsore_wind/data',data_filename);
calmstd = load(f).calmstd;
x1 = calmstd.u;
x2 = calmstd.u_stdv;


h_m = plot(x1,x2, '.', 'Displayname','Measurement data','color','#0071BC');


hold on
% IEC with  Iref = 0.12, Weibull and Lognormal distribution
Iref = 0.12;
filename = sprintf('res/contour_iec_iref_%.2f_%s', Iref, data_filename);
load(filename);
x_combined = [u, NaN, fliplr(u)];          % Combine x values in forward and reverse order
y_combined = [sigma1, NaN,  fliplr(sigma2)];        % Combine y1 and y2 in forward and reverse order
plot(x_combined, y_combined, 'b+', 'Displayname', 'IEC', 'LineWidth', 1);

% IEC with Iref from wind measurement, Weibull and Lognormal distribution
ti = x2./x1;
index = (x1>=14) & (x1<16);
sig_t = x2(index);
ti_t = ti(index);
Iref = mean(ti_t);
filename = sprintf('res/contour_iec_iref_%.2f_%s', Iref, data_filename);
load(filename);
x_combined = [u, NaN, fliplr(u)];          % Combine x values in forward and reverse order
y_combined = [sigma1, NaN,  fliplr(sigma2)];        % Combine y1 and y2 in forward and reverse order
plot(x_combined, y_combined, 'gx', 'Displayname' ,'IEC (data)', 'LineWidth', 1);


% gmm model
filename = sprintf('res/contour_gmm_%s', data_filename);
load(filename);
x_combined = [u, NaN, fliplr(u)];          % Combine x values in forward and reverse order
y_combined = [sigma1, NaN,  fliplr(sigma2)];        % Combine y1 and y2 in forward and reverse order
plot(x_combined, y_combined, 'r+', 'Displayname' ,'GMM', 'LineWidth', 1);


% 

xPmin = [0 0];
xPmax = [40 4];

xlim([xPmin(1), xPmax(1)])
ylim([xPmin(2), xPmax(2)])
xlabel('{\it u} (m/s)')
ylabel('{\it \sigma}_{\it u} (m/s)')
% title('50-year turbulence (45\circ - 135\circ)')
legend('Location','northwest','FontName','Times New Roman')
box on;
set(gca,'FontSize',14,'FontName','Times New Roman')
% 