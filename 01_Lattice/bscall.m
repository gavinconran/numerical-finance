%%% Black-Scholes Call and Put option function

function [C, P] = bscall(S0, K, r, sigma, T)
    d1 = (log(S0/K) + (r + 0.5*sigma^2)*T) / (sigma*sqrt(T));
    d2 = (log(S0/K) + (r - 0.5*sigma^2)*T) / (sigma*sqrt(T));
    C = S0 * normcdf(d1) - K*exp(-r*T)*normcdf(d2);
    
    % Put-Call parity
    P = C - S0 + K*exp(-r*T);
    % return the Call and Put prices for the European option
    bscall = [C, P];
    
    
    
    