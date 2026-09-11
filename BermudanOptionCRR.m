function Bermudan_price=BermudanOptionCRR(F0,K,T,B,sigma,r,q,M)

% Binomial Tree
% The option has a 4-month TTM and exercises at the end of every month

dt = T / M;
% 1-month step = M/4;
Exer_times=[ceil(M/4);ceil(M/2);ceil(3*M/4)]; % % ceil() picks the first discrete step on or just after the exact time fraction

u = exp(sigma*sqrt(dt));
d = 1/u;
p = (1-d) / (u-d); 

% final nodes generation (expiry T)
j = 0:M; % number of up-moves
FT = F0 * (u.^j) .* (d.^(M - j));

% Payoff 
V = max(FT - K, 0);

% Backward Induction with Bermudian condition 
for i=(M-1):-1:0
    
    % Continuation Value 
    V = p * V(2:end) + (1 - p) * V(1:end-1);
    
    % check 
    if ismember(i, Exer_times)
        t = i*dt; % time in years
        
        % price vector of Forward
        j_current = 0:i;
        F_current = F0 * (u.^j_current) .* (d.^(i - j_current));
        
        % computation of the spot price
        S_current = F_current * exp(-(r - q) * (T - t));
        
        % normalized exercise payoff
        Bt = exp(-r * (T - t));
        exer_payoff = max(S_current - K, 0) / Bt;
        
        % bermudian condition 
        V = max(V, exer_payoff);
    end
end

% final conversion 
Bermudan_price = B * V(1);
