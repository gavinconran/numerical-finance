% 1Lattice. Exercise 14
% Let (S t ) t≥0 follow the Black-Scholes model with r = 5%, σ = 20% and S 0 = 100.
% Consider further a European call option with maturity T = 1 year and strike K = 100, and denote
% C 0 (K) its value (to be determined) at inception. Determine a Monte-Carlo procedure to estimate
% the value C 0 (K) of this call option. Use MATLAB to output the following graph:

%%% PART 2: with 1000 paths, show the convergence of the estimator with respect to the number of time
%%% subintervals;


clear all; close all;
graphics_toolkit("gnuplot");
pkg load financial


%%% PART 1
function plot = MonteCarloEstimator2(S0, r, sigma, Tsubs_List, T, sims, K)

% Inputs : S0 - stock price
%        : r - risk free interest rate
%        : sigma - volatility
%        : Tsubs_List - List of number of Time sub intervals
%        : T - time duration
%
% Output: plot - true or false

Price_t_interval = zeros(1, length(Tsubs_List));
index = 1;
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
    Price_t_interval(index) = sum(Est_Price_T)/length(Est_Price_T);
    index += 1;
endfor    

Price_t_interval
% Plot Convergence of Estimation wrt Time intervals   
figure(1)   
plot(Tsubs_List, Price_t_interval);
xlabel('Number of Time Intervals');
ylabel('Stock Price');
str = sprintf('Convergence of Estimation wrt no. of Time Intervals');
title(str); 
ylim([50, 200]);

plot = true;
endfunction  

%%% run part2
S0 = 100;
r=0.05;
sigma=0.20;
T=1.0;
n_min = 1; % min number of time steps
n_max = 10; # 50; 500; % max number of time steps
delta_n = 1; # 5; 50;
% Compute number of steps
steps = ((n_max-n_min)/delta_n)+1;
% Create a vector of number of steps
Tsubs_List = linspace(n_min,n_max, steps); 
sims = 1000;
K = 100;
MonteCarloEstimator2(S0, r, sigma, Tsubs_List, T, sims, K);