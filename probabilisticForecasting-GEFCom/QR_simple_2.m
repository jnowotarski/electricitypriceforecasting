function[QR_hh] = QR_simple_2(Y_hh, X_hh, X_f, Tau_vec)


n=length(Tau_vec);
options= optimset('MaxFunEvals',100000,'MaxIter',50000);

A_initial=inv(X_hh'*X_hh)*X_hh'*Y_hh;
QR_hh = [];

rho_hh = [];
for nn=1:n
    tau = Tau_vec(nn);
    %a_tau = 1-(2*(tau-0.5))^2;
    rho=@(r)sum(r*tau-(r<0).*r);%..*(weights(1:end-1,1));
    A=fminsearch(@(p)rho(Y_hh-X_hh*p),A_initial,options);
    
    rho_help = @(r)sum((r*tau-(r<0).*r));
    
    
    QR_hh = [QR_hh X_f*A];
    
    
    A_initial = A;
end

for hh=1:24
    QR_hh(hh,:) = sort(QR_hh(hh,:));
end
%QR_hh = sort(QR_hh);

end

