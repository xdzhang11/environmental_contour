% Comparison of marginal distributions
% Copulas (nataf) has a predefined marginal theoretical distributions, 
% thus comparison is not needed.
clear;
clc;
close all;
addpath('../utilities');
data_filename = 'offshore_detrend';  % Data file to load
load('res/dis_measurement.mat');
load('res/data_measurement.mat');

% (3)--Gaussian mixture model
sample_title = 'GMM';
% sample_filename = sprintf('res/sample_gmm_%s.mat', data_filename);
% load(sample_filename);
model_filename = sprintf('res/model_gmm_%s.mat', data_filename);
gmModel = load(model_filename).gmModel;
pd = gmm_df(gmModel);
mxpdf = pd.mxpdf;
mxpdf2 = pd.mxpdf2;
pdfx1 = pd.pdfx1;
poex1 = pd.poex1;
pdfx2 = pd.pdfx2;
poex2 = pd.poex2;

xPmin = [0 0];    % lower bound for figure plot
xPmax = [40 4];   % upper bound for figure plot

% fit measurement data to theoretical distributions
% --x1 to Weibull
pd_f1 = fitdist(x1, 'Weibull'); %fitted pd1
pd1_wb = pdf('Weibull',x_pd1,pd_f1.A,pd_f1.B);
% --x2 to lognormal
mu_j = mean(x2);
sigma_j = std(x2);
ln_mu_j = log(mu_j^2/sqrt(sigma_j^2+mu_j^2));
ln_sigma_j = sqrt(log(sigma_j^2/mu_j^2+1));
pd2_ln = pdf('Lognormal',x_pd2,ln_mu_j,ln_sigma_j);


%% Marginal pdf comparison of x1
figure1=figure('Position', [100, 100, 1120, 420]);
    for k = 1:2
        subplot(1,2,k)
        switch k
            case 1
                plot(x_pd1,pd1,'k','DisplayName', 'measurement','LineWidth',1.5);
            case 2
                semilogy(x_pd1, pd1,'k','DisplayName', 'measurement','LineWidth',1.5);
        end
        hold on
        plot(x_pd1,pd1_wb,'b--', 'DisplayName', 'Weibull','LineWidth',1.5)
        plot(x_pd1, pdfx1(x_pd1),'r-.', 'Displayname',sample_title,'LineWidth',1.5)
        xlabel('{\it u} (m/s)')
        ylabel('Probability density')
        legend()
        set(gca,'FontSize',14,'FontName','Times New Roman')
        hold off
    end

%% Marginal poe comparison of x1
figure2 = figure('Position', [100, 100, 560, 420]);
    semilogy(x_icd1, icd1,'k','DisplayName', 'measurement','LineWidth',1.5);
    yticks([1e-6  1e-5 1e-4 1e-3 1e-2 1e-1])
    yticklabels({'10^{-6}','10^{-5}','10^{-4}','10^{-3}','10^{-2}','10^{-1}'})
    hold on
    plot(x_icd1, poex1(x_icd1),'r-.', 'Displayname','GMM','LineWidth',1.5);
    xlabel('{\it u} (m/s)')
    ylabel('Probability of exceedance')
    legend()
    xlim([xPmin(1), xPmax(1)])
    ylim([1e-6,1])
    box on;
    set(gca,'FontSize',16,'FontName','Times New Roman')
    hold off


%% Marginal pdf comparison of x2
figure3=figure('Position', [100, 100, 1120, 420]);
    for k = 1:2
        subplot(1,2,k)
        switch k
            case 1
                plot(x_pd2,pd2,'k','DisplayName', 'measurement','LineWidth',1.5);
            case 2
                semilogy(x_pd2, pd2,'k','DisplayName', 'measurement','LineWidth',1.5);
        end
        hold on
        plot(x_pd2,pd2_ln,'b--', 'DisplayName', 'Lognormal','LineWidth',1.5)
        plot(x_pd2, pdfx2(x_pd2),'r-.', 'Displayname',sample_title,'LineWidth',1.5)
        xlabel('{\it \sigma}_{\it u} (m/s)')
        ylabel('Probability density')
        legend()
        set(gca,'FontSize',14,'FontName','Times New Roman')
        hold off
    end

%% Marginal poe comparison of x2
figure4 = figure('Position', [100, 100, 560, 420]);
    semilogy(x_icd2,icd2,'k','DisplayName', 'Measurement','LineWidth',1.5);
    hold on
    plot(x_icd2, poex2(x_icd2),'r-.', 'Displayname','GMM','LineWidth',1.5)
    xlabel('{\it \sigma}_{\it u} (m/s)')
    ylabel('Probability density')
    legend()
    xlim([xPmin(2), xPmax(2)])
    ylim([1e-6,1])
    box on;
    set(gca,'FontSize',16,'FontName','Times New Roman')
    hold off