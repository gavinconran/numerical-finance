% BinAmerCall.m
% Computes the Binomial tree call option pricing for a particular N

function price = BinAmerCall(S0, K, r, sigma, T, n); 
% Function to calculate the price of a vanilla American
% Call option using a Cox Ross Rubinstein binomial tree.

% Inputs : S0 - stock price
%        : K - strike
%        : r - risk free interest rate
%        : sigma - volatility
%        : T - time duration
%        : n - number of steps
%
% Output: price - the option price

% Calculate the Cox Ross Rubinstein model parameters
dt = T / n; % time step 
a = exp(r*dt);          % discount rate
u = exp(sigma*sqrt(dt));  % up
d = 1/u;                % down
p = (a-d)/(u-d);        % martingale probability


% Compute the Cox Ross Rubinstein underlying price tree
lattice = zeros(n+1,n+1);
lattice(1,1) = S0;
for i = 2:n+1
    lattice(1:i-1,i) = lattice(1:i-1,i-1)*u;
    lattice(i,i) = lattice(i-1,i-1)*d;
end

% Calculate the value at expiry
valueLattice = zeros(size(lattice));
valueLattice(:,end) = max(K-lattice(:,end),0);

% Loop backwards to get values at the earlier times
steps = size(lattice,2)-1;
for idx = steps:-1:1
    valueLattice(1:idx,idx) = ...
        exp(-r*dt)*(p*valueLattice(1:idx,idx+1) ...
        + (1-p)*valueLattice(2:idx+1,idx+1));
        valueLattice(1:idx,idx) = ...
                    max(lattice(1:idx,idx)-K,valueLattice(1:idx,idx));
end

% Output the option price
price = valueLattice(1); #valueTree(1);      

endfunction


%%% TEST
T=1.1
S0=100
sigma = 0.2
K=95
r=0.03
n=5

result = BinAmerCall(S0, K, r, sigma, T, n)

