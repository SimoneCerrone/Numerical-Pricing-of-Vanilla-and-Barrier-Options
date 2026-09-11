function Vega=VegaKO(F0,K,KO,B,T,sigma,N,flagNum)

dsigma=0.01; % Volatility bump

if flagNum==1
    % CRR tree approach - central finite difference method
    P_up = EuropeanOptionKOCRR(F0, K, KO, B, T, sigma + dsigma, N);
    P_down = EuropeanOptionKOCRR(F0, K, KO, B, T, sigma - dsigma, N);
    Vega =dsigma* (P_up - P_down) / (2 * dsigma);
    
elseif flagNum==2
    % MC approach - central finite difference method
    rng(1); % fix the seed, in order not to model noise
    P_up = EuropeanOptionKOMC(F0, K, KO, B, T, sigma + dsigma, N);
    rng(1); 
    P_down = EuropeanOptionKOMC(F0, K, KO, B, T, sigma - dsigma, N);
    Vega= dsigma* (P_up - P_down) / (2 * dsigma);

elseif flagNum==3
    % closed formula
    % By linearity of the Greek, we compute the Vega of the Call options
    % and of the Digital option
    d1_K = (log(F0 / K) + 0.5 * sigma^2 * T) / (sigma * sqrt(T));
    d1_KO = (log(F0 / KO) + 0.5 * sigma^2 * T) / (sigma * sqrt(T));
    d2_KO = d1_KO - sigma * sqrt(T);

    Vega_K=B*F0*sqrt(T)*exp(-(d1_K)^2/2)/sqrt(2*pi)*dsigma;
    Vega_KO=B*F0*sqrt(T)*exp(-(d1_KO)^2/2)/sqrt(2*pi)*dsigma;
    Vega_Digital = - B * normpdf(d2_KO) * (d1_KO / sigma)*dsigma;
    
    % Linearity, according to the European barrier pricing formula
    Vega=Vega_K-Vega_KO-(KO-K)*Vega_Digital;
else 
    error('No valid flagNum. Choose 1 (CRR), 2 (MC) o 3 (Exact).');
end
end