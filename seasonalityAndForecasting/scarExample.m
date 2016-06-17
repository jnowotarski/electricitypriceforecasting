
data = load('NPdata_2013-2016.txt');
startd = 1;
endd = 360*24;
Ndays = 7;

% calculate 7 day-ahead predictions for mSCARX with wavelets, J=10
arxFlag = 2; % multi-day ARX
seasFlag = 1; % wavelets for seasonal decomposition
smoothing = 10; % wavelet decomposition level
example1 = scar(data,Ndays,startd,endd,arxFlag,seasFlag,smoothing);
% calculate WMAE for the prediction:
wmae=100*mean(abs(example1(:,3)-example1(:,4)))/mean(example1(:,3))

% calculate 7 day-ahead predictions for SCARX with wavelets, lambda=1e11
arxFlag = 2; % ARX
seasFlag = 2; % HP filter for seasonal decomposition
smoothing = 1e11; % smoothing parameter lambda in HP filter
example2 = scar(data,Ndays,startd,endd,arxFlag,seasFlag,smoothing);
wmae=100*mean(abs(example2(:,3)-example2(:,4)))/mean(example2(:,3))
