function optionPrice=EuropeanOptionCRR(F0,K,B,T,sigma,N,flag)

% starting parameters for CRR
dt=T/N;
u=exp(sigma*sqrt(dt));
d=1/u;

q=(1-d)/(u-d); % martingale probability under the T-Forward measure

% Final values vector
FN=F0.*(u.^(0:N).*d.^(N:-1:0)); 

% Choose Call or Put
if(flag==1)
    V=max(FN-K,0);
elseif flag==-1
    V=max(K-FN,0);
else
    error('Not valid: use 1 for Call, -1 for Put.');
end

% Bacward induction
for i = N:-1:1
    V = q * V(2:i+1) + (1 - q) * V(1:i);   
end

optionPrice=V*B;
end

