
[<img src="https://github.com/QuantLet/Styleguide-and-FAQ/blob/master/pictures/banner.png" width="880" alt="Visit QuantNet">](http://quantlet.de/index.php?p=info)

## [<img src="https://github.com/QuantLet/Styleguide-and-Validation-procedure/blob/master/pictures/qloqo.png" alt="Visit QuantNet">](http://quantlet.de/) **QRA** [<img src="https://github.com/QuantLet/Styleguide-and-Validation-procedure/blob/master/pictures/QN2.png" width="60" alt="Visit QuantNet 2.0">](http://quantlet.de/d3/ia)

```yaml

Name of Quantlet : QRA

Published in : Computational Statistics

Description : 'Calculates a series of day ahead interval forecasts of electricity spot price. The
QRA method uses a set of individual point forecasts and the output is a pair of quantile forecasts.
Based on: J. Nowotarski, R. Weron (2015) Computing electricity spot price prediction intervals
using quantile regression and forecast averaging, Computational Statistics 30(3), 791-803 (doi:
10.1007/s00180-014-0523-0). '

Keywords : 'electricity price, forecasting, interval forecast, probabilistic forecast, quantile
regression, forecast combination, forecast averaging'

Author : Jakub Nowotarski

Submitted : Wed, June 15 2016 by Jakub Nowotarski

Datafile : NPdata_2013-2016.txt

```


### MATLAB Code:
```matlab
function IntFor = qra(y,X,tau)
%QRA Perform Quantile Regression Averaging
%   QRA(Y,X,TAU) returns an interval forecast at confidence level 
%   (1-tau)*100% obtained using QUANTILE REGRESSION AVERAGING (QRA) [1]. 
%   In the quantile regression model Y is an independent variable in and 
%   variable of interest to have prediction intervals constructed, X is a 
%   matrix of individual forecasts corresponding to values of Y, i.e. real
%   observation and their forecasts for a particular time have the same 
%   indices. 1-TAU is the confidence level. The forecasts is computed for H
%   steps ahead, where H is difference between X's and Y's lengths (has to
%   be positive integer).
%
%   [1] J. Nowotarski, R. Weron (2014) Computing electricity spot price 
%   prediction intervals using quantile regression and forecast averaging,
%   Computational Statistics, forthcoming. Working paper version available 
%   from: http://ideas.repec.org/p/wuu/wpaper/hsc1312.html.
%   Written by Jakub Nowotarski and Rafal Weron (2014.07.30)

% Preliminaries
N = size(X,2);
T = length(y);
if length(X)-T<1
    error('X matrix must contain predictions of future y`s values.')
end

% lower bound
w = quantreg(X(1:T,:),y,tau/2);
if N>1
    IntFor(:,1) = X(T+1:end,:)*w;
else
    IntFor(:,1) = [ones(length(X)-T,1) X(T+1:end,:)]*w;
end

% upper bound
w = quantreg(X(1:T,:),y,1-tau/2);
if N>1
    IntFor(:,2) = X(T+1:end,:)*w;
else
    IntFor(:,2) = [ones(length(X)-T,1) X(T+1:end,:)]*w;
end


function [p,stats]=quantreg(x,y,tau,order,Nboot);
% Quantile Regression
% 
% USAGE: [p,stats]=quantreg(x,y,tau[,order,nboot]);
% 
% INPUTS: 
%   x,y: data that is fitted. (x and y should be columns)
%        Note: that if x is a matrix with several columns then multiple
%        linear regression is used and the "order" argument is not used.
%   tau: quantile used in regression. 
%   order: polynomial order. (default=1)
%          (negative orders are interpreted as zero intercept)
%   nboot: number of bootstrap surrogates used in statistical inference.(default=200)
%
% stats is a structure with the following fields:
%      .pse:    standard error on p. (not independent)
%      .pboot:  the bootstrapped polynomial coefficients.
%      .yfitci: 95% confidence interval on "polyval(p,x)" or "x*p"
%
% [If no output arguments are specified then the code will attempt to make a default test figure
% for convenience, which may not be appropriate for your data (especially if x is not sorted).]
%
% Note: uses bootstrap on residuals for statistical inference. (see help bootstrp)
% check also: http://www.econ.uiuc.edu/~roger/research/intro/rq.pdf
%
% EXAMPLE:
% x=(1:1000)';
% y=randn(size(x)).*(1+x/300)+(x/300).^2;
% [p,stats]=quantreg(x,y,.9,2);
% plot(x,y,x,polyval(p,x),x,stats.yfitci,'k:')
% legend('data','2nd order 90th percentile fit','95% confidence interval','location','best')
% 
% For references on the method check e.g. and refs therein:
% http://www.econ.uiuc.edu/~roger/research/rq/QRJEP.pdf
%
%  Copyright (C) 2008, Aslak Grinsted

%   This software may be used, copied, or redistributed as long as it is not
%   sold and this copyright notice is reproduced on each copy made.  This
%   routine is provided as is without any express or implied warranties
%   whatsoever.


if nargin<3
    error('Not enough input arguments.');
end
if nargin<4, order=[]; end
if nargin<5, Nboot=200; end

if (tau<=0)|(tau>=1),
    error('the percentile (tau) must be between 0 and 1.')
end

if size(x,1)~=size(y,1)
    error('length of x and y must be the same.');
end

if numel(y)~=size(y,1)
    error('y must be a column vector.')
end

if size(x,2)==1
    if isempty(order)
        order=1;
    end
    %Construct Vandermonde matrix.
    if order>0
        x(:,order+1)=1;
    else
        order=abs(order);
    end
    x(:,order)=x(:,1); %flipped LR instead of 
    for ii=order-1:-1:1
        x(:,ii)=x(:,order).*x(:,ii+1);
    end
elseif isempty(order)
    order=1; %used for default plotting
else
    error('Can not use multi-column x and specify an order argument.');
end


pmean=x\y; %Start with an OLS regression

rho=@(r)sum(abs(r-(r<0).*r/tau));

p=fminsearch(@(p)rho(y-x*p),pmean);

if nargout==0
    [xx,six]=sortrows(x(:,order));
    plot(xx,y(six),'.',x(six,order),x(six,:)*p,'r.-')
    legend('data',sprintf('quantreg-fit ptile=%.0f%%',tau*100),'location','best')
    clear p
    return
end 

if nargout>1
    %calculate confidence intervals using bootstrap on residuals
    
    yfit=x*p;
    resid=y-yfit;
    
    stats.pboot=bootstrp(Nboot,@(bootr)fminsearch(@(p)rho(yfit+bootr-x*p),p)', resid);
    stats.pse=std(stats.pboot);
    
    qq=zeros(size(x,1),Nboot);
    for ii=1:Nboot
        qq(:,ii)=x*stats.pboot(ii,:)';
    end
    stats.yfitci=prctile(qq',[2.5 97.5])';
    
end
```
