function result=startols(data,Ndays,startd,endd,alphaci,window)
% Compute combined forecast with OLS method.
% See also J. Nowotarski, E. Raviv, S. Trück, R. Weron (2014) 
%          An empirical comparison of alternate schemes for combining 
%          electricity spot price forecasts, Energy Economics 46, 395-412 
%          (doi: 10.1016/j.eneco.2014.07.014)

Naci = length(alphaci);
result = zeros(Ndays*24,8*Naci+4);
pricecap=1e6;
weights=[];

for j = 1:Ndays
    % display current day and number of days to be forecasted
    disp([j Ndays])
    
% initialize 'result' matrix
    result((j-1)*24+1:j*24,1:3) = data(endd+(j-1)*24+1:endd+j*24,1:3); % result((j-1)*24+hour,1:3) = data(endd+(j-1)*24+hour,1:3);

    % perform one-day-ahead forecast
    [resultaux,w] = forecastols(data(startd:endd+j*24,:),alphaci,window); % data(startd:endd+(j-1)*24+hour,:),alphaci,window);
    result((j-1)*24+1:j*24,4) = min(pricecap, resultaux(:,1)); % result((j-1)*24+hour,4) = min(pricecap, resultaux(:,1));
    result((j-1)*24+1:j*24,5) = resultaux(:,1); % result((j-1)*24+hour,5) = resultaux(:,1);

    % correct results for price cap
    for i=1:Naci
        result((j-1)*24+1:j*24,(i-1)*8+6) = min(pricecap, resultaux(:,(i-1)*4+2)); % result((j-1)*24+hour,(i-1)*8+6) = min(pricecap, resultaux(:,(i-1)*4+2));
        result((j-1)*24+1:j*24,(i-1)*8+7) = min(pricecap, resultaux(:,(i-1)*4+3)); % result((j-1)*24+hour,(i-1)*8+7) = min(pricecap, resultaux(:,(i-1)*4+3));
        result((j-1)*24+1:j*24,(i-1)*8+8) = min(pricecap, resultaux(:,(i-1)*4+4)); % result((j-1)*24+hour,(i-1)*8+8) = min(pricecap, resultaux(:,(i-1)*4+4));
        result((j-1)*24+1:j*24,(i-1)*8+9) = min(pricecap, resultaux(:,(i-1)*4+5)); % result((j-1)*24+hour,(i-1)*8+9) = min(pricecap, resultaux(:,(i-1)*4+5));
        result((j-1)*24+1:j*24,(i-1)*8+(10:13)) = resultaux(:,(i-1)*4+(2:5)); % result((j-1)*24+hour,(i-1)*8+(10:13)) = resultaux(:,(i-1)*4+(2:5));
    end;
    
end

end

function [prediction,weights]=forecastols(data,alphaci,window)

alphaci = (1+alphaci/100)/2;
Naci = length(alphaci);
prediction = zeros(1,4*Naci+1);

for hour=1:24
    
    % preliminaries
    price = data(hour:24:end-24,3);
    ARmodels = data(hour:24:end,4:end);
    if ~strcmp(window,'ext')
        price = price(end+1-window:end);
        ARmodels = ARmodels(end-window:end,:);
    end
    X=ARmodels(1:end-1,:);
    
% perform point forecasts
weights(hour,:) = regress(price,[ones(length(price),1) ARmodels(1:end-1,:)]);
ols = [ones(length(price)+1,1) ARmodels]*weights(hour,:)';
prediction(hour,1) = ols(end);
    
% perform interval forecast
Npred = length(ols)-1;
auxpred = sort(price-ols(1:end-1));
for i=1:Naci
    prediction(hour,(i-1)*4+2)=prediction(hour,1)+(auxpred(ceil(Npred*(1-alphaci(i)))));
    prediction(hour,(i-1)*4+3)=prediction(hour,1)+(auxpred(floor(Npred*alphaci(i))));
    prediction(hour,(i-1)*4+4)=prediction(hour,1)-(norminv(alphaci(i))*std(auxpred));
    prediction(hour,(i-1)*4+5)=prediction(hour,1)+(norminv(alphaci(i))*std(auxpred));
end
end

end