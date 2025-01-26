clear;
clc;
addpath('../utilities');

%% get measurement data
calmstd = load('../data/hovsore_wind/data/offshore_detrend.mat').calmstd;
x1 = calmstd.u;
x2 = calmstd.u_stdv;
save('res/data_measurement.mat', 'x1','x2');


%--------------------------------------------------------------------------
%% marginal distrbution fitting
% marginal distribution of x1
[x_pd1,pd1,x_cd1,cd1,x_icd1,icd1] = epf(x1, 100); % histogram method

    % fit x1 to lognormal distribution 
    mu_i = mean(x1);
    sigma_i = std(x1);
    ln_mu_i = log(mu_i^2/sqrt(sigma_i^2+mu_i^2));
    ln_sigma_i = sqrt(log(sigma_i^2/mu_i^2+1));
    pd1_ln = pdf('Lognormal',x_pd1,ln_mu_i,ln_sigma_i);
    % fit x1 to Weibull distribution
    pd_f1 = fitdist(x1, 'Weibull'); %fitted pd1
    pd1_wb = pdf('Weibull',x_pd1,pd_f1.A,pd_f1.B);
    % [M,V] = wblstat(pd_f1.A,pd_f1.B);
    % fit x2 to Normal distribution
    pd1_normal = pdf('Normal',x_pd1, mu_i,sigma_i);


% marginal distribution of x2
[x_pd2,pd2,x_cd2,cd2,x_icd2,icd2] = epf(x2, 100); % histogram method

    % fit x2 to lognormal distribution
    mu_j = mean(x2);
    sigma_j = std(x2);
    ln_mu_j = log(mu_j^2/sqrt(sigma_j^2+mu_j^2));
    ln_sigma_j = sqrt(log(sigma_j^2/mu_j^2+1));
    pd2_ln = pdf('Lognormal',x_pd2,ln_mu_j,ln_sigma_j);
    
    % fit x2 to Weibull distribution
    pd_f2 = fitdist(x2, 'Weibull');
    pd2_wb = pdf('Weibull',x_pd2,pd_f2.A,pd_f2.B);
    % [M,V] = wblstat(pd_f2.A,pd_f2.B);

    % fit x2 Normal distribution
    pd2_normal = pdf('Normal',x_pd2, mu_j,sigma_j);

 save('res/dis_measurement.mat', 'x_pd1', 'pd1', 'x_cd1', 'cd1', 'x_icd1', 'icd1', ...
     'x_pd2', 'pd2', 'x_cd2', 'cd2', 'x_icd2', 'icd2');

%--------------------------------------------------------------------------
%% plot marginal distributions
xPmin = [0 0];    % lower bound for figure plot
xPmax = [40 4];   % upper bound for figure plot
figure;
    % plot(x_pd1,pd1,'DisplayName', 'x_1');
    semilogy(x_pd1, pd1,'DisplayName', 'x_1');
    hold on;
    plot(x_pd1,pd1_ln,'-.','DisplayName', 'Lognormal','LineWidth',1.5)
    plot(x_pd1,pd1_wb,':', 'DisplayName', 'Weibull','LineWidth',1.5)
    plot(x_pd1,pd1_normal,'k--', 'DisplayName', 'Gaussian','LineWidth',1.5);
    ylim([1e-6 0.12])
    xlabel('{\it u} (m/s)')
    ylabel('Probability density')
    legend()
    xlim([xPmin(1), xPmax(1)])
    % yticks([1e-6  1e-5 1e-4 1e-3 1e-2 1e-1])
    % yticklabels({'10^{-6}','10^{-5}','10^{-4}','10^{-3}','10^{-2}','10^{-1}'})
    % ylim([1e-7, 10])
    set(gca,'FontSize',14,'FontName','Times New Roman')
    hold off

figure;
    % plot(x_pd2,pd2,'DisplayName', 'x_2');
    semilogy(x_pd2, pd2,'DisplayName', 'x_2');
    hold on;
    plot(x_pd2,pd2_ln,'-.','DisplayName', 'Lognormal','LineWidth',1.5);
    plot(x_pd2,pd2_wb,':', 'DisplayName', 'Weibull','LineWidth',1.5)
    plot(x_pd2,pd2_normal,'k--', 'DisplayName', 'Gaussian','LineWidth',1.5);
    ylim([1e-6 1.4])
    xlabel('\sigma_u (m/2)')
    ylabel('Probability density')
    legend()
    xlim([xPmin(2), xPmax(2)])
    % yticks([1e-6  1e-5 1e-4 1e-3 1e-2 1e-1])
    % yticklabels({'10^{-6}','10^{-5}','10^{-4}','10^{-3}','10^{-2}','10^{-1}'})
    % ylim([1e-7, 10])
    set(gca,'FontSize',14,'FontName','Times New Roman')
    hold off
