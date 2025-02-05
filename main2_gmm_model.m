% Train gmm with written and also inline codes
clear;
clc;
addpath('../gmm');

% Define the filename as a variable
data_filename = 'offshore_detrend.mat';  % Data file to load
model_filename = sprintf('model_gmm_%s', data_filename);  % Modify filename for saving model


% Gaussian mixture model using matlab function
%% load data x1 and x2
f = fullfile('../data/hovsore_wind/data',data_filename);
calmstd = load(f).calmstd;
x1 = calmstd.u;
x2 = calmstd.u_stdv;
X(:,1) = x1;
X(:,2) = x2;

%% GMM training
maxComponent = 8; % number of components

% code written
gmModel = gmm_train(X, maxComponent);


% %% MATLAB inline code
% options = statset('Display','final','MaxIter', 1000);
% AIC = zeros(1,maxComponents);
% gm = cell(1,maxComponents);
% for k = 1:maxComponents
%     gm{k} = fitgmdist(X,k,'Options',options, 'Replicates', 3);
%     AIC(k)= gm{k}.AIC;
% end
% [minAIC,numComponents] = min(AIC);
% gmModel = gm{numComponents};
% 
% % Generate new points
% Xp = random(gmModel, 1e8);
% x_i = Xp(:,1);
% x_j = Xp(:,2);

save(fullfile('res', model_filename),'gmModel');