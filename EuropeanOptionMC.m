function optionPrice=EuropeanOptionMC(F0,K,B,T,sigma,N,flag)

Z=randn(N,1); % We generate N numbers with a gaussian distribution

% FN is simulated using the exact analytical solution of the Geometric Brownian Motion
% We use equivalence in law to simulate the terminal value in a single time step
FN=F0.*exp(-0.5*sigma^2*T+sigma*sqrt(T).*Z); 

if flag==1 % Call
    V=max(FN-K,0);
elseif flag==-1 % Put
    V=max(K-FN,0);
else
    error('Not valido: use 1 for Call, -1 for Put.');
end

optionPrice=mean(V)*B; % Compute the mean and discount it

end