function BarrierOptionPriceExact=EuropeanOptionKOClosed(F0,K,KO,B,TTM,sigma)

M=1; % just to set up a number, we won't use it

% We use the formula written in the report

Call_K = EuropeanOptionPrice(F0,K,B,TTM,sigma,1,M,1);
Call_KO = EuropeanOptionPrice(F0,KO,B,TTM,sigma,1,M,1);
d2 = (log(F0 / KO) - 0.5 * sigma^2 * TTM) / (sigma * sqrt(TTM));
Digital = B*normcdf(d2);

BarrierOptionPriceExact = Call_K - Call_KO - (KO-K)*Digital;