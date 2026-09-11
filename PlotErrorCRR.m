function [M,errorCRR]=PlotErrorCRR(F0,K,B,T,sigma)

m=1:10;
M=2.^m; % Choosing M as requested in the assignment

exact_price = EuropeanOptionPrice(F0,K,B,T,sigma,1,1,1);
for i=1:length(M)
    CRR_price = EuropeanOptionPrice(F0,K,B,T,sigma,2,M(i),1);
    errorCRR(i)=abs(exact_price-CRR_price); % absolute error
end

figure(1)
clf

% CRR error
loglog(M, errorCRR, 'o-', 'LineWidth', 1.5)
hold on

% slope -1
guess = errorCRR(1) * (M(1)./M);
loglog(M, guess, '--', 'LineWidth', 1.5)

% 1 Basis Point
xl = xlim;
loglog(xl, [1e-4 1e-4], ':', 'LineWidth', 2)

grid on
grid minor

xlabel('M steps')
ylabel('numerical error')
title('CRR error scaling [log-log]', 'FontWeight', 'bold')

legend('CRR_{error}','1/M','1 Basisi Point', ...
       'Location','southwest', 'FontSize', 14)

