% Compare conditional pdf and icdf between 1) measurement data; 2) copula sample
% data; 3) gmm
clear;
clc;
addpath('../utilities');
data_filename = 'offshore_detrend';  % Data file to load

%% Input options
% -------------------------------------------------------------------------
% % (1)--Nataf transformation
% sample_title = 'Nataf transformation';
% sample_filename = sprintf('res/sample_nataf_%s.mat', data_filename);
% load(sample_filename);
% [x_pdi,pdi,x_cdi,cdi,x_icdi,icdi] = epf(x_i, 100);
% [x_pdj,pdj,x_cdj,cdj,x_icdj,icdj] = epf(x_j, 100);

% % (2)--Gumbel copula
% sample_title = 'Gumbel copula';
% sample_filename = sprintf('res/sample_copula_gumbel_%s.mat', data_filename);
% load(sample_filename);
% %histogram of the sampled data
% [x_pdi,pdi,x_cdi,cdi,x_icdi,icdi] = epf(x_i, 100);
% [x_pdj,pdj,x_cdj,cdj,x_icdj,icdj] = epf(x_j, 100);


% (3)--Gaussian mixture model
sample_title = 'GMM';
model_filename = sprintf('res/model_gmm_%s.mat', data_filename);
gmModel = load(model_filename).gmModel;
pd = gmm_df(gmModel);
mxpdf = pd.mxpdf;
mxpdf2 = pd.mxpdf2;
pdfx1 = pd.pdfx1;
poex1 = pd.poex1;
pdfx2 = pd.pdfx2;
poex2 = pd.poex2;
cpdfx2 = pd.cpdfx2;

% distribution of the measurement data
load('res/data_measurement.mat');

xPmin = [0 0];    % lower bound for figure plot
xPmax = [40 4];   % upper bound for figure plot

%% Conditional pdf
figure1 = figure('Position', [0, 0, 1120, 560]);
bins = 15:2:27;
numInterval = length(bins);
for k = 1:numInterval-1
    subplot(2,3,k)
    % get the conditional sample at each bin
    index = (x1>=bins(k)) & (x1<bins(k+1));
    x_t = x2(index);
    % get the conditional distribution at current bin
    [x_pd_t,pd_t,~,~,~,~] = epf(x_t, 100);
    plot(x_pd_t,pd_t, 'k-', 'Displayname','Measurement','LineWidth',1.5)

    hold on
    if strcmpi(sample_title, 'GMM')
        ut = bins(k)+1;
        x = [repmat(ut, numel(x_pd_t),1),x_pd_t(:)];
        pd_j = cpdfx2(x(:,1),x(:,2));
    else
        index = (x_i>=bins(k)) & (x_i<bins(k+1));
        x_t = x_j(index);
        [x_pd_t,pd_j,~,~,~,~] = epf(x_t, 100);
    end
    plot(x_pd_t,pd_j, '-.', 'Displayname',sample_title,'LineWidth',1.5)

    hold off
    title(sprintf('$u$ = %d m/s',bins(k)+1),'Interpreter','latex')
    xlabel('{\it \sigma}_{\it u} (m/s)')
    ylabel('Probability density')
    box on;
    set(gca,'FontSize',14,'FontName','Times New Roman')
    if bins(k) == 25
        legend()
    end
end

%% Conditional icdf
figure2 = figure('Position', [0, 0, 1120, 560]);
bins = 15:2:27;
numInterval = length(bins);
% figure('Position', [100 100 1200 800])
for k = 1:numInterval-1
    subplot(2,3,k)
    index = (x1>=bins(k)) & (x1<bins(k+1));
    x_t = x2(index);
    [~,~,~,~,x_icd,icd_t] = epf(x_t, 100);

    semilogy(x_icd,icd_t, 'k-', 'Displayname','Measurement','LineWidth',1.5)
    hold on
    if strcmpi(sample_title, 'GMM')
        ut = bins(k)+1;
        pdfc = @(x) mxpdf2(ut,x)./pdfx1(ut);
        cdfc = @(x) integral(pdfc,x, 10);

        pd_j = zeros(1,numel(x_icd));
        for j = 1:numel(x_icd)
            pd_j(j) = cdfc(x_icd(j));
        end
    else
        index = (x_i>=bins(k)) & (x_i<bins(k+1));
        x_t = x_j(index);
        [~,~,~,~,x_icd,pd_j] = epf(x_t, 100);
    end
    semilogy(x_icd,pd_j, 'r-.', 'Displayname',sample_title,'LineWidth',1.5)

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

