
[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="880" alt="Visit QuantNet">](http://quantlet.de/index.php?p=info)

## [<img src="https://github.com/QuantLet/Styleguide-and-Validation-procedure/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **probabilisticForecast** [<img src="https://github.com/QuantLet/Styleguide-and-Validation-procedure/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/d3/ia)

```yaml

Name of Quantlet : probabilisticForecast

Published in : International Journal of Forecasting

Description : 'Computes a day-ahead probabilistic forecast (99 quantiles) of electricity spot price
for the last task of the price track in the Global Energy Forecasting Competition 2014. Based on K.
Maciejowska, J. Nowotarski (2016) A hybrid model for GEFCom2014 probabilistic electricity price
forecasting. International Journal of Forecasting 32(3), 1051-1056.'

Keywords : 'electricity price, forecasting, probabilistic forecast, quantile regression, forecast
combination, global energy forecasting competition'

Author : Jakub Nowotarski

Submitted : Wed, July 17 2016 by Jakub Nowotarski

Datafile : GEFcom-Task15.txt, Benchmark15.csv, solution15_P.csv

```


### MATLAB Code:
```matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ======== STEP 0: Do some preprocessing ========
clear all; clc

%% Load the data
Y = load('GEFcom-Task15.txt');
Benchmark =  csvread('Benchmark15.csv');
realPrice = csvread('solution15_P.csv',1,2); 
realPrice = realPrice(:,1); % For some versions of matlab

%% Define variables
loadTotal = reshape(Y(:,1),24,size(Y,1)/24)';
loadZonal = reshape(Y(:,2),24,size(Y,1)/24)';
Load = [loadTotal loadZonal];
Price = reshape(Y(:,3),24,size(Y,1)/24)';

% Dummy variable for weekends
T = size(loadTotal,1);
Day = (0:1:T-1)';
Day = Day-floor(Day/7)*7; %sat - 0, sun - 1, mon - 2, ..., fri - 6
Weekend = (Day<2);
Day = [Day  Weekend];
maxLag = 8;

% Dummy variable for peak hours
Peak = zeros(size(Y,1),24);	%Peak hours between 10-19
Peak(:,10:19) = ones(size(Y,1),10);

% For regression models
Y = Price(1:end,:);
Y_forecast = Price(end,:);

% Vector of quantiles
Tau_vec = (1:99)/100;
n = length(Tau_vec);

% Load ratio (for Step 3 - QR)
loadRatio = (Load(maxLag+1:end,1:24)./Load(maxLag:end-1,1:24)).* ...
    (Load(1:end-maxLag,1:24)./Load(2:end-maxLag+1,1:24));

%% Simulations - preliminaries
% Measure the distance of the load from the load pattern of the forecasted day
Load_1 = Load(:,1:24);
Load_1 = mean((Load_1 - repmat(Load_1(end,:),size(Load_1,1),1)).^2,2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ======== STEP 1: Compute point forecast ========
%% Preliminaries
Y_mean_1 = zeros(size(Y,1)-maxLag,24);
Y_mean_2 = zeros(size(Y,1)-maxLag,24);

Load_1 = Load_1(maxLag+1:end,:);

%% Predict 24 hours separately
for hh=1:24
    if hh>5
X_help = [ones(size(Day,1)-maxLag,1) Day(maxLag+1:end,2) log([Y(maxLag:end-1,hh)...
    Y(2:end-maxLag+1,hh) Load(maxLag+1:end,hh)...
    min(Load(maxLag+1:end,1:24),[],2)./max(Load(maxLag+1:end,1:24),[],2)...
    Load(maxLag:end-1,hh) Load(2:end-maxLag+1,hh)  ]) ];
    else
        % For hours 1:4 add P_{T,[23,24]} as explanatory variables
X_help = [ones(size(Day,1)-maxLag,1) Day(maxLag+1:end,2) log([Y(maxLag:end-1,23:24)...
    Y(maxLag:end-1,hh) Y(2:end-maxLag+1,hh)  Load(maxLag+1:end,1:24)...
    min(Load(maxLag+1:end,1:24),[],2)./max(Load(maxLag+1:end,1:24),[],2)...
    Load(maxLag:end-1,hh) Load(2:end-maxLag+1,hh)  ]) ];
    end
    
    Y_help = log(Y(maxLag+1:end,hh));

    % \hat(Y1) = X(1:T+1,:) * \hat(beta)
    Y_mean_1(:,hh) = X_help * (X_help(1:end-1,:)\Y_help(1:end-1,:));
    
    % ARX2: use 10% of days with the most simmilar load patters
    indc = find(Load_1<=quantile(Load_1,0.1));
    % Remove the collumn describing the weekends
    X_help= [X_help(:,1:2) X_help(:,3:end)];	
    X_help_2 = X_help(indc,:);
    Y_help_2 = Y_help(indc,:);
    
    % \hat(Y2) = X(1:T+1,:) * \hat(beta)
    Y_mean_2(:,hh) = X_help * (X_help(1:end-1,:)\Y_help(1:end-1,:));
        
end

%% Return to nominal values
Y_mean_1 = exp(Y_mean_1);
Y_mean_2 = exp(Y_mean_2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ======== STEP 2: Pre-filter the residuals  ========

%% Define variables
T_help = size(Y_mean_1,1);
Y_t = Y(end-T_help+1:end,:);
Day_t = Day(end-T_help+1:end,:);
Load_t = Load(end-T_help+1:end,:);
Load_1_t = Load_1(end-T_help+1:end,:);
Load_2_t = repmat(mean(Load_t')',1,24);	% Average daily load


%% Day type filtering
indc = find(Day_t(:,1)>2); %Pick the day

Y_t = Y_t(indc,:);
Day_t = Day_t(indc,:);
Load_t = Load_t(indc,:);

Load_1_t = Load_1_t(indc,:);
Y_mean_1 = Y_mean_1(indc,:);
Y_mean_2 = Y_mean_2(indc,:);
loadRatio = loadRatio(indc,:);
Peak = Peak(indc,:);

%% Similar load profile filtering
indc = find(Load_1_t<=quantile(Load_1_t,0.1));

Y_t = Y_t(indc,:);
Day_t = Day_t(indc,:);
Load_t = Load_t(indc,:);
Load_1_t = Load_1_t(indc,:);
Y_mean_1 = Y_mean_1(indc,:);
Y_mean_2 = Y_mean_2(indc,:);
loadRatio = loadRatio(indc,:);
Peak = Peak(indc,:);

%% Combine 2 point forecasts
Y_mean_0 = (Y_mean_1 + Y_mean_2)/2;

%% Calculate the residuals
Y_mean_all = repmat(mean(Y_mean_0')',1,24);	% Dialy mean of the prices
Y_hh = Y_t(1:end-1,:)-Y_mean_0(1:end-1,:);	% Residuals

%% Transform matrices into the vectors
Y = reshape(Y_hh',[],1);
loadRatio = reshape(loadRatio',[],1);
Y_mean_1 = reshape(Y_mean_1',[],1);
Y_mean_2 = reshape(Y_mean_2',[],1);
Y_mean_0 = reshape(Y_mean_0',[],1);
Y_mean_all = reshape(Y_mean_all',[],1);
Peak = reshape(Peak',[],1);
Load = reshape(Load_t(:,1:24)',[],1);
Load_2 = reshape(Load_2_t(:,1:24)',[],1);

Y_mean = [Y_mean_1 Y_mean_2];


%% Expected bias filtering
% Normalize Y
Normal_Y = (Y-mean(Y))/std(Y);
Normal_Y_t = (Y(end-23:end,:)-mean(Y))/std(Y);
% a = 4, b = 0
indc = find(and(Normal_Y<0, Normal_Y>-4));




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ======== STEP 3: Quantile regression  ========

%% Preliminiaries
Y = Y(indc,1);
loadRatio = [loadRatio(indc,1); loadRatio(end-23:end,:)];
Y_mean = [Y_mean(indc,:); Y_mean(end-23:end,:)];
Y_mean_0 = [Y_mean_0(indc,:); Y_mean_0(end-23:end,:)];
Y_mean_all = [Y_mean_all(indc,:); Y_mean_all(end-23:end,:)];
Peak = [Peak(indc,1); Peak(end-23:end,:)];
Load = [Load(indc,1); Load(end-23:end,:)];
Load_2 = [Load_2(indc,1); Load_2(end-23:end,:)];

% Y_mean_all = smooth(Y_mean_all,24); % This line is not described in the paper but for some weeks we did use it

T_help = length(Y_mean);
Y_hh = Y;

%% X matrix for quantile regression
X_help = [ones(T_help,1) loadRatio Load Load_2  Y_mean_all  Y_mean(:,1) Y_mean(:,1).^2 ];
X_hh = X_help(1:end-24,:);
X_f = X_help(end-23:end,:);

% Add dummy for peak hours
ind_help = (Peak>0);
X_help_1 = X_help.*repmat(ind_help,1,size(X_help,2));
X_hh = [X_hh X_help_1(1:end-24,:)];
X_f = [X_f X_help_1(end-23:end,:)];

%% Run quantile regression (sorting quantiles is icluded in QR_simple_2)
QR_tt = QR_simple_2(Y_hh, X_hh, X_f, Tau_vec)+ repmat(Y_mean_0(end-23:end,1),1,n);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ======== STEP 4: Post-processing ========

%% Smooth the quantile curves
QR_smoothed = [];
for nn=1:n
    QR_nn = smooth(QR_tt(:,nn),3);
    QR_smoothed = [QR_smoothed QR_nn];
end

%% Plot the probabilistic forecast
% http://www.mathworks.com/help/finance/fanplot.html
historical = [0:24; [Price(end-1,end) realPrice']]';
forecast = [1:24; QR_smoothed']';
fanplot(historical,forecast)
title('Probabilistic forecast')
ylabel('Price')
xlabel('Hours')
set(gca,'xtick',[1:6:24 24])
```
