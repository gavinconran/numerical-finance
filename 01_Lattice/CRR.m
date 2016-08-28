% Convergence of the CRR price to Black-Scholes as a function of the number of nodes N
% using Cox, Ross and Rubenstein (CRR) binomial options pricing model 

clear all; close all;
graphics_toolkit("gnuplot");
pkg load financial

% set   parameters
S0 = 100;
r=0.03;
sigma=0.2;
T=1.1;

%%% k =95
K=95;
n_min = 100; % min number of time steps
n_max = 5000; % max number of time steps
delta_n = 20;

% blsprice is a function of the financial package
%[BS_C,BS_P] = blsprice(S0, K, r, T, sigma);
% compute analytically the price of a European Call and Put option for K=95
[BS_C, BS_P] = bscall(S0, K, r, sigma, T);

% Compute number of steps
steps = ((n_max-n_min)/delta_n)+1;
% Create a vector of number of steps
num_steps = linspace(n_min,n_max, steps); 
% Initialise the vector of Option Prices
BinomPrices = zeros(1, length(num_steps));


% For each number of steps compute the Option Price
index = 1;
for n = num_steps
    BinomPrices(index) = BinEuroCall(S0, K, r, sigma, T, n); 
    %BinomPrices(index) = binPriceCRR(K,S0,r,sigma,delta_T,n,'CALL',false);
    index = index + 1;
end;
 
% Plot for K = 95 
figure(1)                                              
plot(num_steps, BinomPrices,'-', num_steps, ones(1,length(num_steps))*BS_C, 'r-');
xlabel('Number of Nodes (N)');
ylabel('Option Price');
title('Convergence of CRR price to Black-Scholes as a function of number of nodes N');h = legend('Approx Price  ', 'Analytic Price  ');

% k = 100
K=100;
n_min = 100; % min number of time steps
n_max = 3000; % max number of time steps
delta_n = 20;

% blsprice is a function of the financial package
%[BS_C,BS_P] = blsprice(S0, K, r, T, sigma);
% compute analytically the price of a European Call and Put option for K=100
[BS_C, BS_P] = bscall(S0, K, r, sigma, T);

% Compute number of steps
steps = ((n_max-n_min)/delta_n)+1;
% Create a vector of number of steps
num_steps = linspace(n_min,n_max, steps); 
% Initialise the vector of Option Prices
BinomPrices = zeros(1, length(num_steps));

% For each number of steps compute the Option Price
% Compute for even and odd step numbers
index = 1;
for n_even = num_steps
    n_odd = n_even + 1;
    delta_T_even = T / n_even;
    delta_T_odd = T / n_odd;
    BinomPrices_even(index) = BinEuroCall(S0, K, r, sigma, T, n_even); 
    %BinomPrices_even(index) = binPriceCRR(K,S0,r,sigma,delta_T_even, n_even,'CALL',false);
    BinomPrices_odd(index) = BinEuroCall(S0, K, r, sigma, T, n_odd); 
    %BinomPrices_odd(index) = binPriceCRR(K,S0,r,sigma,delta_T_odd, n_odd,'CALL',false);
    index = index + 1;
end;

% Plot for K = 100    
figure(2)                                              
plot(num_steps, BinomPrices_even,'-', num_steps, BinomPrices_odd, '-', num_steps, ones(1,length(num_steps))*BS_C, '-');
xlabel('Number of Nodes (N)');
ylabel('Option Price');
title('Convergence of CRR price to Black-Scholes as a function of number of nodes N');
legend('Approx Price (Even)  ', 'Approx Price (Odd)  ', 'Analytic Price  ');

