clear;
clc
data_filename = 'offshore_detrend';  % Data file to load
f = fullfile('../data/hovsore_wind/data',data_filename);
calmstd = load(f).calmstd;
x1 = calmstd.u;
x2 = calmstd.u_stdv;
X(:,1) = x1;
X(:,2) = x2;

model_filename = sprintf('res/model_gmm_%s.mat', data_filename);
gmModel = load(model_filename).gmModel;
pd = gmm_df(gmModel);
mxpdf = pd.mxpdf;

%% Plot setting
numP = 2500;
xPmin = [0 0];
xPmax = [40 4];
% xPmax = [25 3];
x1Plot = linspace(xPmin(1), xPmax(1), numP);
x2Plot = linspace(xPmin(2), xPmax(2), numP);

%% Pdf plot
figure; 
[xp1,xp2] = meshgrid(linspace(xPmin(1),xPmax(1),sqrt(numP))',linspace(xPmin(2),xPmax(2),sqrt(numP))');
xPlot = [xp1(:) xp2(:)];
y_pdf_t = mxpdf(xPlot);  % Theoretical pdf value
y_pdf_t_r = reshape(y_pdf_t,sqrt(numP),sqrt(numP));   % Reshaped theoretical pdf value
surf(xp1,xp2,y_pdf_t_r);
xlim([xPmin(1), xPmax(1)])
ylim([xPmin(2), xPmax(2)])
xlabel('X_1');
ylabel('X_2');
zlabel('Probability density')
title('Mixture model')


%% Contour plot
figure;
% Contour plot from data
ctL = [1/(365*24*6) 1e-4 1e-3 1e-2]; %contour line values
% Contour plot from histogram
h = histogram2(x1,x2,25,'Normalization','pdf','facecolor','flat');
% Use the middle points of bin edges from histrogram
x1Edges = h.XBinEdges;
xbTemp = x1Edges(1:end-1);
xeTemp = x1Edges(2:end);
x1T = (xbTemp+xeTemp)/2;
x2Edges = h.YBinEdges;
xbTemp = x2Edges(1:end-1);
xeTemp = x2Edges(2:end);
x2T = (xbTemp+xeTemp)/2;
% Get the meshgrid data
[xp1T,xp2T] = meshgrid(x1T,x2T); 
zz = h.Values;  %Z is the density value from histogram plot
zz = zz';        
% Contour Plot
[C,s] = contour(xp1T,xp2T,zz,ctL,'LineStyle','-','color','b','LineWidth',1);
clabel(C,s,'FontSize',12)

% Contour plot from GMM model
hold on 
[xp1,xp2] = meshgrid(linspace(xPmin(1),xPmax(1),sqrt(numP))',linspace(xPmin(2),xPmax(2),sqrt(numP))');
xPlot = [xp1(:) xp2(:)];
%pdf values for contour plot
y_pdf_t = mxpdf(xPlot);  % Theoretical pdf value
y_pdf_t_r = reshape(y_pdf_t,sqrt(numP),sqrt(numP));   % Reshaped theoretical pdf value
% Contour plot from GMM
contour(xp1,xp2,y_pdf_t_r,ctL,'ShowText','on','LineStyle',':','color','r','LineWidth',3);
xlabel('{\it u} (m/s)')
ylabel('\sigma_u (m/s)')
set(gca,'FontSize',14,'FontName','Times New Roman')

