load indFor.mat

alpha = 0.1; % tau=1-alpha is confidence level, forecast interval: [tau/2,1-tau/2]
Ndays = 7;
startd = 1; % beginnng of the calibration period
endd = 182*24; % end of the calibration period

for d=1:Ndays
    for h=1:24
        yh = log(data((d-1)*24+startd-1+h:24:(d-1)*24+endd,3));
        Xh = log(data((d-1)*24+startd-1+h:24:d*24+endd,4:end));

        % QRA - [lower bound, upper bound]
        PI((d-1)*24+h,:) = exp(qra(yh,Xh,alpha));
    end
end

% real price during forecast period
y = data(endd+1:endd+Ndays*24,3);

% calculate coverage
I = PI(:,1)<y & y<PI(:,2);
disp(mean(I))