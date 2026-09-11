function [M,stdEstim]=PlotErrorMC(F0,K,B,T,sigma)

m=1:20;
M=2.^m; % Choosing M as requested in the assignment

for i=1:length(M)
    rng(1);
    Z=randn(M(i),1);
    FN=F0.*exp(-0.5*sigma^2*T+sigma*sqrt(T).*Z);
    V=max(FN-K,0);
    price=B*V;
    stdEstim(i)=std(price)/sqrt(M(i)); % unbiased standard deviation
end

figure(2)
clf

% MCC error
loglog(M, stdEstim, 'o-', 'LineWidth', 1.5)
hold on

% slope -1/2 
guess = stdEstim(1) * (sqrt(M(1))./sqrt(M));
loglog(M, guess, '--', 'LineWidth', 1.5)

% 1 Basis Point
xl = xlim;
loglog(xl, [1e-4 1e-4], ':', 'LineWidth', 2)

grid on
grid minor

xlabel('M simulations')
ylabel('std error')
title('MC error scaling [log-log]', 'FontWeight', 'bold')

legend('MC_{error}','1/\surdM','1 Basis Point', ...
       'Location','southwest','FontSize', 14)


hold off