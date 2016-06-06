function result=startirmse(data,Ndays,startd,endd,alphaci,window)
% Compute combined forecast with IRMSE method.
% See also J. Nowotarski, E. Raviv, S. Trück, R. Weron (2014) 
%          An empirical comparison of alternate schemes for combining 
%          electricity spot price forecasts, Energy Economics 46, 395-412 
%          (doi: 10.1016/j.eneco.2014.07.014)

Naci = length(alphaci);
result = zeros(Ndays*24,8*Naci+4);
pricecap=1e6;

for j = 1:Ndays
    % display current day and number of days to be forecasted
    disp([j Ndays])
    
    % initialize 'result' matrix
    for hour=1:24
        result((j-1)*24+hour,1:3) = data(endd+(j-1)*24+hour,1:3);
    
        % perform one-day-ahead forecast
        resultaux = forecastirmse(data(startd:endd+(j-1)*24+hour,:),alphaci,window);
        result((j-1)*24+hour,4) = min(pricecap, resultaux(:,1));
        result((j-1)*24+hour,5) = resultaux(:,1);

        % correct results for price cap
        for i=1:Naci
            result((j-1)*24+hour,(i-1)*8+6) = min(pricecap, resultaux(:,(i-1)*4+2));
            result((j-1)*24+hour,(i-1)*8+7) = min(pricecap, resultaux(:,(i-1)*4+3));
            result((j-1)*24+hour,(i-1)*8+8) = min(pricecap, resultaux(:,(i-1)*4+4));
            result((j-1)*24+hour,(i-1)*8+9) = min(pricecap, resultaux(:,(i-1)*4+5));
            result((j-1)*24+hour,(i-1)*8+(10:13)) = resultaux(:,(i-1)*4+(2:5));
        end;
    
    end
end



end

function prediction=forecastirmse(data,alphaci,window)
% opis

alphaci = (1+alphaci/100)/2;
Naci = length(alphaci);
prediction = zeros(1,4*Naci+1);

% preliminaries
if strcmp(window,'ext')
    price = data(1:end-1,3);
    ARmodels = data(1:end,4:end);
else
    price = data(end-24*window:end-1,3);
    ARmodels = data(end-24*window:end,4:end);
end
Nmodels = size(data,2)-3;
repprice = repmat(price,1,Nmodels);
    
% perform point forecasts
rmse = sqrt(mean((ARmodels(1:end-1,:)-repprice(1:end,:)).^2));
weights = (1./rmse')/sum(1./rmse);
irmse = ARmodels*weights;
prediction(1) = irmse(end);

% perform interval forecast
Npred = length(irmse)-1;
auxpred = sort(price-irmse(1:end-1));
for i=1:Naci
    prediction((i-1)*4+2)=prediction(1)+(auxpred(ceil(Npred*(1-alphaci(i)))));
    prediction((i-1)*4+3)=prediction(1)+(auxpred(floor(Npred*alphaci(i))));
    prediction((i-1)*4+4)=prediction(1)-(norminv(alphaci(i))*std(auxpred));
    prediction((i-1)*4+5)=prediction(1)+(norminv(alphaci(i))*std(auxpred));
end

end