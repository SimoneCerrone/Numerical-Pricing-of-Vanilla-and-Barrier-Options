% Group 12, AA 2025-2026
%
% Simone Cerrone    (Matr: 303432)
% Mohamad Ternanni  (Matr: 298496)
% Filippo Tacconi   (Matr: 314998)
% Vittoria Tomasini (Matr: 304183)
%
%% Pricing parameters
% All parameters should be put here, in the script and passed to the
% fuctions of interest (generally in a struct)

S0=1;          % underlying price of equity stock  
K=1.1;         % strike price 
r=0.025;       % risk free rate 
TTM=1/3;       % Time to maturity 
sigma=0.212;   % volatility 
d=0.02;        % dividend yield 
nc = 1e6;      % numer of underlying contracts  
Notional = nc*S0;  % Notational 
flag=1;        % flag:  1 call, -1 put

%% Quantity of interest

B=exp(-r*TTM); % Discount factor 

%% a) Pricing the option 

F0=S0*exp(-d*TTM)/B; % Forward in G&C Model
M=100; % M = simulations for MC, steps for CRR;

pricingMode = [1 2 3]; % 1 ClosedFormula, 2 CRR, 3 Monte Carlo
OptionPrice =[];       % vector initialization
for i=1: length(pricingMode)
    rng(1);
    OptionPrice(i) = EuropeanOptionPrice(F0,K,B,TTM,sigma,pricingMode(i),M,flag);
end

exact_price = OptionPrice(1)
CRR_price = OptionPrice(2)
MC_price = OptionPrice(3)
% by increasing M the two prices converge to to the price obtained using the closed formula

fprintf('The total price of the option using closed formula is: %.4f\n',exact_price*Notional);
fprintf('The total price of the option using CRR is: %.4f\n',CRR_price*Notional);
fprintf('The total price of the option using MC formula is: %.4f\n',MC_price*Notional);

%% b) Select M according to the criteria mentioned in the class.

precision=0.0001;  % 1 bp

% CRR
err=200;          % value initialization to enter in the cicle
M_CRR=1;          % statring value of M for CRR
possible=100000;  % value initialization 
counter=0; % counter implemented to check if, once the correct M is found the following values are still lower than 1 bp

% The idea is that since CRR is a discrete time model and the steps can
% either surpass the strike or not depending on the single jump, it might
% happen (as we discover in the next point) that for M time steps the
% error is smaller than 1BP but for M' > M it might turn back to be higher,
% so we wanto to find a value s.t. for the next values all of them remain
% under the treshold (specifically for the next 50 time steps).
% While the error is still above the precision mark the cycle goes on to compute the
% errors: if the error is lower then 1 bp you increase the counter of "good steps" 
% and fix the possible outcome of the algorithm. But if the error goes back up, the counter goes
% back to 0 (since the found result is not a stable) and the possible outcome is resetted. 
% At the end of the cicle, the output found is s.t. for the next 50 iterations it won't 
% come surpass the precision trashold (it is stable).

while (err>= precision || counter<50) 
    CRR_price = EuropeanOptionPrice(F0,K,B,TTM,sigma,2,M_CRR,1);
    err = abs(exact_price - CRR_price);
    if err < precision              % the number of time steps works
        counter = counter + 1;   % update the counter
        possible = min(possible,M_CRR); % this is a possible value for the request
    else              % the error goeas back to being higher than 1BP, the result is unstable
        counter = 0;      % reset the counter
        possible = 10000; % reset the possible value
    end
    M_CRR = M_CRR + 1; % increasing of the number of time steps linearly
end
M_CRR = possible   % possible outcome that we found is stable

% MC
M_MC=0;      % starting value of M for MC
err=200;        % value initialization to enter in the cicle

while(err >= precision)
    M_MC=M_MC + 1000;   % increase M by 1000, otherwise algorithm NOT computationally feasible 
    rng(1);
    Z=randn(M_MC,1);
    FN=F0.*exp(-0.5*sigma^2*TTM+sigma*sqrt(TTM).*Z);
    V=max(FN-K,0);
    price=B*V;
    err=std(price)/sqrt(M_MC);
end % The error is estimated as the unbiased standard deviation of the MC estimator
M_MC

% Using the values of M found, we get:
format long   % to see the actual difference
exact_price = EuropeanOptionPrice(F0,K,B,TTM,sigma,1,M,flag)
CRR_price = EuropeanOptionPrice(F0,K,B,TTM,sigma,2,M_CRR,flag)
MC_price = EuropeanOptionPrice(F0,K,B,TTM,sigma,3,M_MC,flag)
format short

%% c) Errors rescaling with M approximately as 1/𝑀 for CRR and as 1/√𝑀 for MC.

[nCRR,errCRR] = PlotErrorCRR(F0,K,B,TTM,sigma);
M2_CRR = nCRR(find(errCRR < precision, 5))  % values of M_CRR following the powers of 2 convention 
% M=32 accidentally hits the tolerance, but M=64 fails

rng(1) % We fix the seed in order to reproduce always the same graph
[nMC,stdEstim] = PlotErrorMC(F0,K,B,TTM,sigma);
M2_MC = nMC(find(stdEstim < precision, 1)) % value of M_MC following the powers of 2 convention 

%% d) European Barrier Option 

KO = 1.4; % knock-out value 

% CRR pricing 
BarrierOptionPriceCRR = EuropeanOptionKOCRR(F0,K, KO,B,TTM,sigma,M_CRR);

fprintf('Price of the European Barrier Option using CRR: %f \n',BarrierOptionPriceCRR);

% MC pricing 
rng(1);
BarrierOptionPriceMC = EuropeanOptionKOMC(F0,K, KO,B,TTM,sigma,M_MC);

fprintf('Price of the European Barrier Option using MC: %f \n',BarrierOptionPriceMC );

% Closed formula pricing: obtained considering 2 Call options (with different strikes) and a Digital option
Call_K = EuropeanOptionPrice(F0,K,B,TTM,sigma,1,M,1);
Call_KO = EuropeanOptionPrice(F0,KO,B,TTM,sigma,1,M,1);
d2 = (log(F0 / KO) - 0.5 * sigma^2 * TTM) / (sigma * sqrt(TTM));
Digital = B*normcdf(d2);
BarrierOptionPriceExact = Call_K - Call_KO - (KO-K)*Digital;

fprintf('Price of the European Barrier Option using closed formula: %f \n',BarrierOptionPriceExact);

%% e) Plot Vega

% By linearity of the greek and using the closed formula found above:

S0nodes = linspace(0.65,1.45,100); % discrtization of price range 
for i=1:length(S0nodes)
    F0nodes = S0nodes(i)*exp(-d*TTM)/B;
    Vega_CRR(i) = VegaKO(F0nodes,K,KO,B,TTM,sigma,M_CRR,1); % Vega using CRR
    rng(1)
    Vega_MC(i) = VegaKO(F0nodes,K,KO,B,TTM,sigma,M_MC,2);  % Vega using MC
    Vega_exact(i) = VegaKO(F0nodes,K,KO,B,TTM,sigma,1,3);  % Vega using closed formula 
end
figure(3);
plot(S0nodes,Vega_exact,'LineWidth', 1.5);
hold on
plot(S0nodes,Vega_CRR,'LineWidth', 1.5);
plot(S0nodes,Vega_MC,'LineWidth', 1.5);
grid on 
grid minor
xlabel('Underlying price')
ylabel('\nu','FontSize',14)
title('\nu of European Up-and-Out Call','FontWeight', 'bold')
legend('Closed Formula','CRR','MC','Location','southwest','FontSize', 14)
hold off

%% f)  American Barrier

% Price
S0nodes = linspace(0.65,1.45,100); % discrtization of price range
for i=1:length(S0nodes)

    % European Barrier - closed formula
    F0nodes = S0nodes(i)*exp(-d*TTM)/B;
    Call_K = EuropeanOptionPrice(F0nodes,K,B,TTM,sigma,1,M,1);
    Call_KO = EuropeanOptionPrice(F0nodes,KO,B,TTM,sigma,1,M,1);
    d2 = (log(F0nodes / KO) - 0.5 * sigma^2 * TTM) / (sigma * sqrt(TTM));
    Digital = B*normcdf(d2);
    EUR_barr(i) = Call_K-Call_KO-(KO-K)*Digital;
    
    % American Barrier - closed formula
    AM_barr(i)=EuropeanOptionAmericanKOClosed(S0nodes(i),K, KO,B,TTM,sigma,r,d);
end
figure(4)
plot(S0nodes,EUR_barr,'LineWidth', 1.5);
hold on
plot(S0nodes,AM_barr,'LineWidth', 1.5);
xline(KO, '--');
grid on 
grid minor
xlabel('Underlying price')
ylabel('Options Price')
title('Price Comparison','FontWeight', 'bold')
legend('EU Barrier Option ','AM Barrier Option','Barrier value','Location','northwest','FontSize', 14)
hold off

% For the original parameters
AM_barr_exact=EuropeanOptionAmericanKOClosed(S0,K, KO,B,TTM,sigma,r,d);
fprintf('\nPrice of the European barrier option: %.4f\n',BarrierOptionPriceExact);
fprintf('Price of the American barrier option: %.4f\n',AM_barr_exact);

% Delta
dS=0.0001; % Price bump
for i=1:length(S0nodes)

    % European Barrier - closed formula
    F0nodes = S0nodes(i)*exp(-d*TTM)/B;
    d1_K = (log(F0nodes / K) + 0.5 * sigma^2 * TTM) / (sigma * sqrt(TTM));
    d1_KO = (log(F0nodes / KO) + 0.5 * sigma^2 * TTM) / (sigma * sqrt(TTM));
    d2_KO = d1_KO - sigma * sqrt(TTM);
    Delta_K = exp(-d * TTM) * normcdf(d1_K);
    Delta_KO = exp(-d * TTM) * normcdf(d1_KO);
    Delta_Digital = B * normpdf(d2_KO) / (S0nodes(i) * sigma * sqrt(TTM));
    Delta_EUR(i) = Delta_K - Delta_KO - (KO - K) * Delta_Digital;
    
    % American Barrier - central finite difference method
    F0_up = (S0nodes(i) + dS) * exp(-d * TTM) / B;
    F0_down = (S0nodes(i) - dS) * exp(-d * TTM) / B;
    P_up_exact=EuropeanOptionAmericanKOClosed(S0nodes(i)+dS,K, KO,B,TTM,sigma,r,d);
    P_down_exact=EuropeanOptionAmericanKOClosed(S0nodes(i)-dS,K, KO,B,TTM,sigma,r,d);
    Delta_AM(i) = (P_up_exact - P_down_exact) / (2 * dS);
end
figure(5)
plot(S0nodes,Delta_EUR,'LineWidth', 1.5);
hold on
plot(S0nodes,Delta_AM,'LineWidth', 1.5);
xline(KO, '--');
grid on 
grid minor
xlabel('Underlying price')
ylabel('\Delta','FontSize',14)
title('\Delta Comparison','FontWeight', 'bold')
legend('EU Barrier Option ','AM Barrier Option','Barrier value','Location','northwest','FontSize', 14)
hold off

% Vega
dsigma = 0.01; % Sigma bump
for i=1:length(S0nodes)

    % European Barrier - closed formula
    F0nodes = S0nodes(i)*exp(-d*TTM)/B;
    Vega_EUR(i) = VegaKO(F0nodes,K,KO,B,TTM,sigma,M,3);
    
    % American Barrier - central finite difference method
    P_AM_vol_up = EuropeanOptionAmericanKOClosed(S0nodes(i), K, KO, B, TTM, sigma + dsigma, r,d);
    P_AM_vol_down = EuropeanOptionAmericanKOClosed(S0nodes(i), K, KO, B, TTM, sigma - dsigma, r,d);
    Vega_AM(i) = (P_AM_vol_up - P_AM_vol_down) / 2;
end
figure(6)
plot(S0nodes,Vega_EUR,'LineWidth', 1.5);
hold on
plot(S0nodes,Vega_AM,"LineWidth",1.5);
xline(KO, '--');
grid on 
grid minor
xlabel('Underlying price')
ylabel('\nu','FontSize',14)
title('\nu Comparison','FontWeight', 'bold')
legend('EU Barrier Option ','AM Barrier Option','Barrier value','Location','southwest','FontSize', 14)
hold off

%% g) Antithetic Variables techinique 

[nMC,stdEstim1] = PlotErrorMC(F0,K,B,TTM,sigma);
hold on
[M_anti,stdEstimANTI] = PlotErrorANTIMC(F0,K,B,TTM,sigma);
grid on
legend('stdEstim','1/sqrt(M)','1BP','stdEstimANTI','Location','southwest');
title('MC error scaling with antithetic variables tecnique[log-log]', 'FontWeight', 'bold')

% Error comparison with M found in point b)

% Standard MC
rng(1);
Z=randn(M_MC,1);
FN=F0.*exp(-0.5*sigma^2*TTM+sigma*sqrt(TTM).*Z);
V=max(FN-K,0);
price_MC=B*V;
stdEstim_MC=std(price_MC)/sqrt(M_MC)

% Antithetic MC
rng(1)
Z=randn(floor(M_MC/2),1);

FN_plus=F0.*exp(-0.5*sigma^2*TTM+sigma*sqrt(TTM).*Z);
V_plus=max(FN_plus-K,0);
price_plus=B*V_plus;

FN_minus=F0.*exp(-0.5*sigma^2*TTM+sigma*sqrt(TTM).*(-Z));
V_minus=max(FN_minus-K,0);
price_minus=B*V_minus;
    
price_pairs=(price_plus + price_minus) / 2;
stdEstim_ANTI=std(price_pairs)/sqrt(floor(M_MC/2))
% Commented in the associated function (PlotErrorANTIMC)

% The error decreases even by simulating half of the random numbers
% We reduced the computational cost and improved the performance

%% h) Bermudan Option

Berm_price = BermudanOptionCRR(F0,K,TTM,B,sigma,r,d,M)

%% i) Dividend yield comparison 

div=linspace(0,0.05,100); % dividend yield
for i=1:length(div)
    F0_new = S0*exp(-div(i)*TTM)/B; % 
    Bermudan(i) = BermudanOptionCRR(F0_new,K,TTM,B,sigma,r,div(i),M_CRR);
    European(i) = EuropeanOptionPrice(F0_new,K,B,TTM,sigma,2,M_CRR,1); %  Price using CRR for consistency
end
figure(7)
plot(div,European,'LineWidth',1.5);
hold on
plot(div,Bermudan,'LineWidth',1.5);
grid on 
grid minor
xlabel('dividend yield')
ylabel('Options price')
legend('European Option','Bermudan Option','Location','southwest','FontSize', 14);