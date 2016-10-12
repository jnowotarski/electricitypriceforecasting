% Perform forecast combination for electricity prices.
% Algorithm uses extending calibration window for weights estimation.
% See also J. Nowotarski, E. Raviv, S. Trück, R. Weron (2014) 
%          An empirical comparison of alternate schemes for combining 
%          electricity spot price forecasts, Energy Economics 46, 395-412 
%          (doi: 10.1016/j.eneco.2014.07.014)

load indResults.mat

startd = 1; % starting day for weights estimation
endd = startd+28*24-1; % last day for weights estimation
Ndays = 7; % numer of days to be forecasted


% SIMPLE AVERAGING
y = data(:,3);
X = data(:,4:end);
pred(:,1) = mean(X(endd+1:endd+24*Ndays,:),2);

% Models with weight estimation
% IRMSE
result=startirmse(data,Ndays,startd,endd,[],'ext');
pred(:,2) = result(:,4);
% OLS
result=startols(data,Ndays,startd,endd,[],'ext');
pred(:,3) = result(:,4);
% CLS
result=startcls(data,Ndays,startd,endd,[],'ext');
pred(:,4) = result(:,4);
