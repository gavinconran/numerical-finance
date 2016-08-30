% 1Lattice. Exercise 14
% Let (S t ) t≥0 follow the Black-Scholes model with r = 5%, σ = 20% and S 0 = 100.
% Consider further a European call option with maturity T = 1 year and strike K = 100, and denote
% C 0 (K) its value (to be determined) at inception. Determine a Monte-Carlo procedure to estimate
% the value C 0 (K) of this call option. Use MATLAB to output the following graph:

%%% PART 3: with 1000 paths and 100 time subintervals, plot the function K → C 0 (K) for K = 20, . . . , 150.

clear all; close all;
graphics_toolkit("gnuplot");
pkg load financial

function plot = MonteCarloEstimator(S0, r, sigma, Tsubs_List, T, sims, K_List)

% Inputs : S0 - stock price
%        : r - risk free interest rate
%        : sigma - volatility
%        : Tsubs_List - List of number of Time sub intervals
%        : T - time duration
%
% Output: plot - true or false
for steps = Tsubs_List
    % Create a vector of number of steps
    time_steps = linspace(0, 1, steps+1); 
    % Initialise the vector of Estimated Stock Price at time T
    Price_t = zeros(sims, length(time_steps));
    Price_t(:,1) = S0;
    t = T / (steps+1);
    
    % Initialise the vector of Estimated Stoke Prices
    Est_Price_T = zeros(1, sims);
    num_sims = linspace(1,sims, sims);
    
    Price_sum = 0; #S0;
    for i=1:sims #num_sims     
        % For each time step compute the Stock Price St 
        for j = 2:length(time_steps)         
            Price_t(i,j) = Price_t(i,j-1)*exp((r-sigma.^2/2)*t + sigma*sqrt(t)*normrnd(0,1));
        endfor
        % Compute Estimator for i simulations
        Price_sum += Price_t(i,end);
        Est_Price_T(i) = Price_sum / i;
    endfor 
endfor    

CallOptionPrice_K = zeros(1, length(K_List));
index = 1;
for K=K_List
    CallOptionPrice_K(index++) = max(Est_Price_T(end)-K,0);
end    
    
% Plot C(K) as  afunction of K   
figure(1)   
plot(K_List, CallOptionPrice_K);
xlabel('Strike Price (K)');
ylabel('Call Option Price');
str = sprintf('Call Option Price C(K) in relation to Strike Price (K)');
title(str); 
ylim([-10, 110]);

plot = true;
endfunction  
  
     


S0 = 100;
r=0.05;
sigma=0.20;
T=1.0;
Tsubs_List = [100];
sims = 1000;

% compute K_List
K_min = 20; 
K_max = 150; 
delta_K = 10;
steps = ((K_max-K_min)/delta_K) + 1;
% Create a vector of number of steps
K_List = linspace(K_min,K_max, steps); 
%run simulation 
MonteCarloEstimator(S0, r, sigma, Tsubs_List, T, sims, K_List, 3);


