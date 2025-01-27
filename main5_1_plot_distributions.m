% 1) data points comparison
% 2&3) Compare marginal distributions from 1)measurement data; 2)theoretical
% distributions; 3)sampled data
% 4&5) Compare conditional pdf and icdf between 1)measurement data: 2) sample
% data
clear;
clc;
addpath('../utilities');

%% Input options
% -------------------------------------------------------------------------
% (1)--Nataf transformation
% sample_title = 'Nataf transformation';
% load('data/sample_nataf.mat');

% (2)--Gumbel copula
sample_title = 'Gumbel copula';
load('res/sample_copula_gumbel.mat');


% distribution of the measurement data
load('res/dis_measurement.mat');
load('res/data_measurement.mat');

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

% histogram of the sampled data
[x_pdi,pdi,x_cdi,cdi,x_icdi,icdi] = epf(x_i, 100);
[x_pdj,pdj,x_cdj,cdj,x_icdj,icdj] = epf(x_j, 100);


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


%% marginal distribution comparison of x1
figure2=figure('Position', [100, 100, 1120, 420]);
for k = 1:2
    subplot(1,2,k)
    switch k
        case 1
            plot(x_pd1,pd1,'DisplayName', 'measurement','LineWidth',1.5);
        case 2
            semilogy(x_pd1, pd1,'DisplayName', 'measurement','LineWidth',1.5);
    end
    hold on
    plot(x_pd1,pd1_wb,'--', 'DisplayName', 'Weibull','LineWidth',1.5)
    plot(x_pdi, pdi,'+', 'Displayname',sample_title,'LineWidth',1.5)
    xlabel('x_1')
    ylabel('Probability density')
    legend()
    set(gca,'FontSize',14,'FontName','Times New Roman')
    hold off
end


%% Marginal distribution comparison of x2
figure3=figure('Position', [100, 100, 1120, 420]);
for k = 1:2
    subplot(1,2,k)
    switch k
        case 1
            plot(x_pd2,pd2,'DisplayName', 'measurement','LineWidth',1.5);
        case 2
            semilogy(x_pd2, pd2,'DisplayName', 'measurement','LineWidth',1.5);
    end
    hold on
    plot(x_pd2,pd2_ln,'--', 'DisplayName', 'Lognormal','LineWidth',1.5)
    plot(x_pdj, pdj,'+', 'Displayname',sample_title,'LineWidth',1.5)
    xlabel('x2')
    ylabel('Probability density')
    legend()
    set(gca,'FontSize',14,'FontName','Times New Roman')
    hold off
end


%% Conditional pdf
figure4 = figure('Position', [0, 0, 1680, 840]);
bins = 3:2:27;
numInterval = length(bins);
for k = 1:numInterval-1
    subplot(3,4,k)
    % get the conditional sample at each bin
    index = (x1>=bins(k)) & (x1<bins(k+1));
    x_t = x2(index);
    % get the conditional distribution at current bin
    [x_pd_t,pd_t,~,~,~,~] = epf(x_t, 100);
    plot(x_pd_t,pd_t, '-', 'Displayname','Measurement','LineWidth',1.5)

    hold on
    index = (x_i>=bins(k)) & (x_i<bins(k+1));
    x_t = x_j(index);
    [x_pd_t,pd_t,~,~,~,~] = epf(x_t, 100);
    plot(x_pd_t,pd_t, '-.', 'Displayname',sample_title,'LineWidth',1.5)    
    
    hold off
    title(sprintf('$u$ = %d m/s',2*k+2),'Interpreter','latex')
    xlabel('{\it \sigma}_{\it u} (m/s)')
    ylabel('Probability density')
    box on;
    set(gca,'FontSize',14,'FontName','Times New Roman')
    if k ==10
    legend()
    end
end

%% Conditional icdf
figure5 = figure('Position', [0, 0, 1120, 560]);
bins = 15:2:27;
numInterval = length(bins);
% figure('Position', [100 100 1200 800])
for k = 1:numInterval-1
    subplot(2,3,k)
    index = (x1>=bins(k)) & (x1<bins(k+1));
%     if sum(index)<100
%         continue
%     end
    x_t = x2(index);
    [~,~,~,~,x_icd_t,icd_t] = epf(x_t, 100);

    semilogy(x_icd_t,icd_t, '-', 'Displayname','Measurement','LineWidth',1.5)
    hold on

    index = (x_i>=bins(k)) & (x_i<bins(k+1));
%     if sum(index)<100
%         continue
%     end
    x_t = x_j(index);
    [~,~,~,~,x_icd,icd] = epf(x_t, 100);
    semilogy(x_icd,icd, '-.', 'Displayname',sample_title,'LineWidth',1.5)
    
    hold off
    title(sprintf('$u$ = %d m/s',bins(k)+1),'Interpreter','latex')
    xlabel('{\it \sigma}_{\it u} (m/s)')
    ylabel('Probability of exceedance')
    yticks([1e-5 1e-4 1e-3 1e-2 1e-1 1])
    yticklabels({'10^{-5}','10^{-4}','10^{-3}','10^{-2}','10^{-1}','1'})
    box on;
    set(gca,'FontSize',14,'FontName','Times New Roman')
    if bins(k) ==25
    legend()
    end
end

