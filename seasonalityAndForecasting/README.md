
[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="880" alt="Visit QuantNet">](http://quantlet.de/index.php?p=info)

## [<img src="https://github.com/QuantLet/Styleguide-and-Validation-procedure/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **seasonalityAndForecasting** [<img src="https://github.com/QuantLet/Styleguide-and-Validation-procedure/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/d3/ia)

```yaml

Name of Quantlet : seasonalityAndForecasting

Published in : Energy Economics

Description : 'Computes an example of forecast for electricity spot prices with the (m)SCARX model.
(m)SCARX calculates a day-ahead prediction using decoposition of the raw electricity prices into a
seasonal component and remaining stochastic residuals. Based on J. Nowotarski, R. Weron (2016) On
the importance of the long-term seasonal component in day-ahead electricity price forecasting,
Energy Economics 57, 228-235 (doi: 10.1016/j.eneco.2016.05.009).'

Keywords : electricity price, forecasting, seasonality, wavelets, hodrick-prescot filter, ARX

Author : Jakub Nowotarski

Submitted : Wed, June 15 2016 by Jakub Nowotarski

Datafile : NPdata_2013-2016.txt

```


### MATLAB Code:
```matlab
% Calculate day-ahead predictions for the 7 days following 
% the 360-day calibration window (GEFCom2014 dataset) 
% Requires:
%   scar.m - SCAR-type model estimation and forecasting routine
%   

data = load('NPdata_2013-2016.txt');
startd = 1; % first hour of the calibration window
endd = 360*24; % last hour of the calibration window
Ndays = 104*7; % number of days to be predicted

% mSCARX-S_10
arxFlag = 2; % mSCARX
seasFlag = 1; % wavelets for seasonal decomposition
smoothing = 10; % wavelet decomposition level
example1 = scar(data,Ndays,startd,endd,arxFlag,seasFlag,smoothing);

% mSCARX-HP_5e8
arxFlag = 2; % mSCARX
seasFlag = 2; % HP filter for seasonal decomposition
smoothing = 5e8; % smoothing parameter lambda in HP filter
example2 = scar(data,Ndays,startd,endd,arxFlag,seasFlag,smoothing);

% SCARX-S_9
arxFlag = 1; % SCARX
seasFlag = 1; % wavelets for seasonal decomposition
smoothing = 9; % wavelet decomposition level
example3 = scar(data,Ndays,startd,endd,arxFlag,seasFlag,smoothing);

% SCARX-HP_1e8
arxFlag = 1; % SCARX
seasFlag = 2; % HP filter for seasonal decomposition
smoothing = 1e8; % smoothing parameter lambda in HP filter
example4 = scar(data,Ndays,startd,endd,arxFlag,seasFlag,smoothing);

% calculate mean WMAE:
for ii=168:168:Ndays*24
    TT = ii-167:ii;
    wmae1(ii/168) = 100*mean(abs(example1(TT,3)-example1(TT,4)))/mean(example1(TT,3));
    wmae2(ii/168) = 100*mean(abs(example2(TT,3)-example2(TT,4)))/mean(example2(TT,3));
    wmae3(ii/168) = 100*mean(abs(example3(TT,3)-example3(TT,4)))/mean(example3(TT,3));
    wmae4(ii/168) = 100*mean(abs(example4(TT,3)-example4(TT,4)))/mean(example4(TT,3));
end    
    
disp(['                WMAE'])
disp(['mSCARX-S_10     ' num2str(mean(wmae1))])
disp(['mSCARX-HP_5e8   ' num2str(mean(wmae2))])
disp(['SCARX-S_9       ' num2str(mean(wmae3))])
disp(['SCARX-HP_1e8    ' num2str(mean(wmae4))])


```
