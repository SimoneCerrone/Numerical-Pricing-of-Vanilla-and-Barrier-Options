function optionPrice=EuropeanOptionAmericanKOClosed(S0,K, KO,B,TTM,sigma,r,d)

% if the spot price is above the barrier, the option is null and void
if (S0 >= KO)
    optionPrice = 0;
    return;
end

% otherwise:
F0 = S0 * exp(-d*TTM) / B; % Forrward price
F1=(KO^2/S0)*exp(-d*TTM)/B; % Forward price based on the reflected spot (KO^2 / S0)

% chopped payoff method
option1=EuropeanOptionKOClosed(F0,K,KO,B,TTM,sigma);
option2=EuropeanOptionKOClosed(F1,K,KO,B,TTM,sigma);

optionPrice = option1 - (KO/S0)^(2*(r-d-0.5*sigma^2)/sigma^2)*option2;
end