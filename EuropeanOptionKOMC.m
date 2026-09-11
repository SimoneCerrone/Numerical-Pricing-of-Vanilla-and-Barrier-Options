function optionPrice=EuropeanOptionKOMC(F0,K, KO,B,T,sigma,N)

% same pricing method as the european call option
Z=randn(N,1);
FN=F0.*exp(-0.5*sigma^2*T+sigma*sqrt(T).*Z);

V=max(FN-K,0);

% we check if FN reaches barrier
V(FN >= KO) = 0;

optionPrice=mean(V)*B;

end