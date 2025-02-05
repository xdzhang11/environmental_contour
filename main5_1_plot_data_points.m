clear;
clc;
addpath('../utilities');

data_filename = 'offshore_detrend.mat';  % Data file to load

% % (1)--Nataf transformation
% sample_title = 'Nataf transformation';
% sample_filename = sprintf('res/sample_nataf_%s.mat', data_filename);
% load(sample_filename);

% % (2)--Gumbel copula
% sample_title = 'Gumbel copula';
% sample_filename = sprintf('res/sample_copula_gumbel_%s.mat', data_filename);
% load(sample_filename);

% (3)--Gaussian mixture model
sample_title = 'GMM';
sample_filename = sprintf('res/sample_gmm_%s.mat', data_filename);
load(sample_filename);


xPmin = [0 0];    % lower bound for figure plot
xPmax = [40 4];   % upper bound for figure plot

% distribution of the measurement data
load('res/data_measurement.mat');


%% data points plot
figure1=figure('Position', [100, 100, 1120, 420]);
subplot(1,2,1);
scatter(x1, x2,10,'.');
xlim([xPmin(1), xPmax(1)])
ylim([xPmin(2), xPmax(2)])
xlabel('{\it u} (m/s)')
ylabel('{\it \sigma}_{\it u} (m/s)')
title('Measurement')
box on;
set(gca,'FontSize',14,'FontName','Times New Roman')

subplot(1,2,2)
scatter(x_i,x_j,10,'c.') % Scatter plot with points of size 10
xlim([xPmin(1), xPmax(1)])
ylim([xPmin(2), xPmax(2)])
xlabel('{\it u} (m/s)')
ylabel('{\it \sigma}_{\it u} (m/s)')
title(sample_title)
box on;
set(gca,'FontSize',14,'FontName','Times New Roman')