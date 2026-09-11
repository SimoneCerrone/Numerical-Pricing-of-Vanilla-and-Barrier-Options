function [M,stdEstim]=PlotErrorANTIMC(F0,K,B,T,sigma)

% It mirrors the structure of PlotErrorMC but implements the 
% variance reduction technique

m=1:20;
M=2.^m; % Choosing M as requested in the assignment

for i=1:length(M)
    rng(1) % Fix the seed for reproducibility of the convergence plot
    Z=randn(floor(M(i)/2),1); % Generate half of the required random numbers

    % By symmetry:

    FN_plus=F0.*exp(-0.5*sigma^2*T+sigma*sqrt(T).*Z);
    V_plus=max(FN_plus-K,0);
    price_plus=B*V_plus;

    FN_minus=F0.*exp(-0.5*sigma^2*T+sigma*sqrt(T).*(-Z));
    V_minus=max(FN_minus-K,0);
    price_minus=B*V_minus;
    
    price_pairs=(price_plus + price_minus) / 2; % Average the paired payoffs
    % Compute standard error
    stdEstim(i)=std(price_pairs)/sqrt(floor(M(i)/2));
end
loglog(M,stdEstim,'LineWidth',1.5);
end