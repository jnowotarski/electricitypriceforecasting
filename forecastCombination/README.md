
[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="880" alt="Visit QuantNet">](http://quantlet.de/index.php?p=info)

## [<img src="https://github.com/QuantLet/Styleguide-and-Validation-procedure/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **forecastCombination** [<img src="https://github.com/QuantLet/Styleguide-and-Validation-procedure/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/d3/ia)

```yaml

Name of Quantlet : forecastCombination

Published in : Energy Economics

Description : 'Computes an example of forecast combination for electricity spot prices. Uses 4
methods: simple average, OLS averaging restricted regression and IRMSE averaging. Based on J.
Nowotarski, E. Raviv, S. Trueck, R. Weron (2014) An empirical comparison of alternate schemes for
combining electricity spot price forecasts, Energy Economics 46, 395-412 (doi:
10.1016/j.eneco.2014.07.014).'

Keywords : electricity price, forecasting, forecast combination, forecast averaging

Author : Jakub Nowotarski

Submitted : Mon, June 6 2016 by Jakub Nowotarski

Datafile : indResults.mat

```


### MATLAB Code:
```matlab
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

```
