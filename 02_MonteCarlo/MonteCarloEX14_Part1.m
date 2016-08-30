% 1Lattice. Exercise 14
% Let (S t ) t≥0 follow the Black-Scholes model with r = 5%, σ = 20% and S 0 = 100.
% Consider further a European call option with maturity T = 1 year and strike K = 100, and denote
% C 0 (K) its value (to be determined) at inception. Determine a Monte-Carlo procedure to estimate
% the value C 0 (K) of this call option. Use MATLAB to output the following graph:

%%% PART 1: with 10 time subintervals, show the convergence of the estimator 
%%% with respect to the number of simulated paths;


clear all; close all;
graphics_toolkit("gnuplot");
pkg load financial

function plot = MonteCarloEstimator(S0, r, sigma, Tsubs_List, T, sims, K)

% Inputs : S0 - stock price
%        : r - risk free interest rate
%        : sigma - volatility
%        : Tsubs_List - List of number of Time sub intervals
%        : T - time duration
%        : sims - number of simulations
%        : K - Strike Price
%
% Output: plot - true or false

% compute simulations
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

%Plot Price at time T histogram
figure(1)
nbins=25;
hist(Price_t(:,end), nbins);
xlabel('Stock Price at time T');
ylabel('Number of Stock');
str = sprintf('Price at time T histogram');
title(str); 

% Plot Convergence of Estimation wrt no. Simulated Paths   
figure(2)   
plot(num_sims, Est_Price_T);
xlabel('Number of Simulations (Sims)');
ylabel('Stock Price');
str = sprintf('Convergence of Estimation wrt no. Simulated Paths');
title(str); 
ylim([50, 200]);

plot = true;
endfunction  

%%% run part1
S0 = 100;
r=0.05;
sigma=0.20;
T=1.0;
Tsubs_List = [10];
sims = 10000;
K=100;
MonteCarloEstimator(S0, r, sigma, Tsubs_List, T, sims, K);