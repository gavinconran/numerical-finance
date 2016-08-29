% 1Lattice. Exercise 14
% Let (S t ) t≥0 follow the Black-Scholes model with r = 5%, σ = 20% and S 0 = 100.
% Consider further a European call option with maturity T = 1 year and strike K = 100, and denote
% C 0 (K) its value (to be determined) at inception. Determine a Monte-Carlo procedure to estimate
% the value C 0 (K) of this call option. Use MATLAB to output the following graph:
%%% PART 1: with 10 time subintervals, show the convergence of the estimator with respect to the number
%%%  of simulated paths;
%%% PART 2: with 1000 paths, show the convergence of the estimator with respect to the number of time
%%% subintervals;
%%% PART 3: with 1000 paths and 100 time subintervals, plot the function K → C 0 (K) for K = 20, . . . , 150.

clear all; close all;
graphics_toolkit("gnuplot");
pkg load financial


%%% PART 1
function plot = MonteCarloEstimator(S0, r, sigma, Tsubs_List, T, sims)

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
    t = T / steps;
    
    % Initialise the vector of Estimated Stoke Prices
    Est_Price_T = zeros(1, sims);
    num_sims = linspace(1,sims, sims);
    
    Price_sum = 0;
    for i=num_sims     
        % For each time step compute the Stock Price St 
        for j = 2:length(time_steps)            
            Price_t(i,j) = Price_t(i,j-1)*exp((r-sigma^2/2)*t + sigma*sqrt(t)*normrnd(0,1));
        endfor
        %Price_t
        Price_sum += Price_t(i,end);
        Est_Price_T(i) = Price_sum / i;
        %Price_sum / i
    endfor 
endfor    

% Plot    
figure(1)   
plot(num_sims, Est_Price_T);
xlabel('Number of Simulations (Sims)');
ylabel('Stock Price');
str = sprintf('Convergence of Estimation wrt no. Simulated Paths');
title(str); 
                                      
Est_Price_T(end)
plot = true;

endfunction  

% PART 1
S0 = 100;
r=0.05;
sigma=0.20;
T=1.0;
Tsubs_List = [10]; #[10, 100, 1000];
sims = 1000;
K=100;
MonteCarloEstimator(S0, r, sigma, Tsubs_List, T, sims);


