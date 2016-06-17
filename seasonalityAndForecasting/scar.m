function result = scar(data,Ndays,startd,endd,arxFlag,seasFlag,smoothing)
% SCAR: Calculates day-ahead predictions of electricity spot price given
% historical prices and forecasted loads with the (m)SCAR(X) model.
%
%   Input data:
%       DATA - 4-column matrix (date, hour, price, load forecast),
%       NDAYS - number of days to be forecasted,
%       STARTD - index of the first day of the calibration period,
%       ENDD - index of the last day (in the first step) of the calibration 
%           period,
%       ARXFLAG - a flag for ARX model selection:
%           ARXFLAG=1 --> ARX
%           ARXFLAG=2 --> multi-day ARX
%       SEASFLAG - a flag for deseasonalize method selection:
%           SEASFLAG=0 --> use raw data for computing the prediction (no
%               seasonality removal)
%           SEASFLAG=1 --> wavelets db24
%           SEASFLAG=2 --> HP filter
%       SMOOTHING - wavelet decomposition level (if SEASFLAG==1) or lambda
%           (smoothing parameter) in HP filter (if SEASFLAG==2)
%   Output format: 4-column matrix (date, hour, price, forecasted price). 
%
%   Reference(s): 
%   [1] J. Nowotarski, R. Weron (2016) On the importance of the long-term 
%       seasonal component in day-ahead electricity price forecasting, 
%       submitted.

result = zeros(Ndays*24,4);

for j = 1:Ndays
    % Uncomment to display current day and number of days to be predicted
    % disp([j Ndays])
       
    % initialize 'result' matrix
    result((j-1)*24+1:j*24,1:3) = data(endd+(j-1)*24+1:endd+j*24,1:3);
    
    % perform one-day-ahead forecast
    resultaux = forecastscar(data(startd+(j-1)*24:endd+j*24,:),arxFlag,seasFlag,smoothing);
    result((j-1)*24+1:j*24,4) = resultaux(:,1);
    
end;


function prediction=forecastscar(data,arxFlag,seasFlag,smoothing)
%FORECASTSCAR Internally used by SCAR.
%   PREDICTION=FORECASTSCAR(...) returns one-day-ahead point forecasts. 
%   PREDICTION is a vector with 24 length.
%
%   Input data:
%       DATA - 4-column matrix (date, hour, price, load forecast),
%       ARXFLAG - a flag for ARX model selection:
%           ARXFLAG=1 --> ARX
%           ARXFLAG=1 --> multi-day ARX
%       SEASFLAG - a flag for deseasonalize method selection:
%           SEASFLAG=0 --> use raw data for computing the prediction (no
%           seasonality removal)
%           SEASFLAG=1 --> wavelets db24
%           SEASFLAG=2 --> HP filter
%       SMOOTHING - wavelet decomposition level (if SEASFLAG==1) or lambda
%       (smoothing parameter) in HP filter (if SEASFLAG==2)

% weekday of starting point
i = weekday(datenum(num2str(data(1,1)),'yyyymmdd'))-1;

N = length(data);
data = [data zeros(N,3)];
for j=1:24:N
    switch mod(i,7)
        case 6
            data(j:j+23,7)=1;   % Saturday
        case 0
            data(j:j+23,8)=1;   % Sunday
        case 1
            data(j:j+23,9)=1;   % Monday
    end;
    i=i+1;
end;

% initialize 'prediction' matrix
prediction = zeros(24,1);

% Deseasonalize
if seasFlag==1 % wavelets
    dwtmode('sp0','nodisp')
    ptemp = [data(1:end-24,3);mean(data(end-47:end-24,3))]; % so that sp0 extension is on the correct level
    [C,L] = wavedec(log(ptemp),smoothing+2,'db24');
    ltsc = wrcoef('a',C,L,'db24',smoothing); ltsc = ltsc(1:end-1,:);
elseif seasFlag==2 % HP filter
  	ltsc = hpfilter(log(data(1:end-24,3)),smoothing);    
end
    data(:,5) = data(:,3); % real price
    data(:,6) = [ltsc;ltsc(end-23:end,:)]; % LTSC with prediction
    data(:,3) = log(data(:,3)) - data(:,6); % residuals

% compute forecasts
for hour = 1:24
    % preliminaries
    price = data(hour+168:24:end-24,3);
    price_168 = data(hour:24:end-168,3);
    price_min = data(169-24:end-24,3);
    price_min = reshape(price_min,24,length(price_min)/24);
    price_min = min(price_min)';
    
    if ~seasFlag % no seasonality in the model - remove the mean
        price = log(price);
        mc = mean(price);
        price = price-mc;
        price_168 = log(price_168);
        price_168 = price_168-mean(price_168);
        price_min = log(price_min);
        price_min = price_min-mean(price_min);
    else
        ltsc = data(hour+168:24:end-24,6);
    end

    loadr = data(hour+168:24:end,4);
    loadr = log(loadr);
        
    % calibrate ARX/mARX model
    D = [data(hour+168:24:end,7:9) ones(length(data(hour+168:24:end,7:9)),1)];
    y = price;
    
    if arxFlag==1 % ARX
        X = [price(3:end-1)...
            price(2:end-2) price_168(4:end-1)...
            price_min(4:end-1) loadr(4:end-1) D(4:end-1,1:3)];
        X_fut = [price(end)...
            price(end-1) price_168(end)...
            price_min(end) loadr(end) D(end,1:3)];
    elseif arxFlag==2 % multi-day ARX
        X = [repmat(price(3:end-1),1,4).*D(4:end-1,:)...
            price(1:end-3).*D(4:end-1,3)...
            price(2:end-2) price_168(4:end-1)...
            price_min(4:end-1) loadr(4:end-1) D(4:end-1,:)];
        X_fut = [price(end)*D(end,:)...
            price(end-2)*D(end,3)...
            price(end-1) price_168(end)...
            price_min(end) loadr(end) D(end,:)];
    end
    
    beta = regress(y(4:end),X);
    prog = X_fut*beta;
    
    if seasFlag
        prediction(hour,1) = exp(prog+ltsc(end));
    else
        prediction(hour,1) = exp(prog+mc);
    end
    
end