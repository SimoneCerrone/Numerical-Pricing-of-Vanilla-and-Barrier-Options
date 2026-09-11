function optionPrice=EuropeanOptionKOCRR(F0,K, KO,B,T,sigma,N)

% same pricing method as the european call option
dt=T/N;
u=exp(sigma*sqrt(dt));
d=1/u;

q=(1-d)/(u-d);

% Final values vector
FN=F0.*(u.^(0:N).*d.^(N:-1:0)); 

V=max(FN-K,0);

% we check if FN reaches barrier
V(FN >= KO) = 0;

% Bacward induction
for i = N:-1:1
    V = q * V(2:i+1) + (1 - q) * V(1:i);
end

optionPrice=V*B;

end