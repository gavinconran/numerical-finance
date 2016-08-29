% 1Lattice. Exercise 5 
% Exercise 5. Consider a European call option with maturity T > 0 and strike K > 0 evaluated
% at time zero, written on a stock price process (S t ) t≥0 following the Black-Scholes model. 
% Consider the following values: S 0 = 100, r = 0 and σ = 25%. Plot the convergence
% of the CRR tree when K < 90, K = 100 and K = 110. Plot also the convergence of the tree for
% even and odd increasing number of time steps. Comment the obtained results.

clear all; close all;
graphics_toolkit("gnuplot");
pkg load financial

function plot = latticeEx1(S0, r, sigma, K_List, n_min, n_max, delta_n, T)

% Inputs : S0 - stock price
%        : r - risk free interest rate
%        : sigma - volatility
%        : K_List - List of Strike prices
%        : n_min - minimum number of steps used in simulation
%        : n_max - mamimum number of steps used in simulation
%        : delta_n - step size
%        : T - time duration
%
% Output: plot - true or false



%%% k < 90
%K=85;
fig_num = 1; % update figure number for plotting
for K = K_List
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

    % Plot for K = 90    
    figure(fig_num)                                              
    plot(num_steps, BinomPrices_even,'-', num_steps, BinomPrices_odd, '-', num_steps, ones(1,length(num_steps))*BS_C, '-');
    xlabel('Number of Nodes (N)');
    ylabel('Option Price');
    str = sprintf('Convergence of CRR price to BSM as a function of N. K = %d', K);
    title(str); 
    %title('Convergence of CRR price to BSM as a function of N. K < %d', K);
    legend('Approx Price (Even)  ', 'Approx Price (Odd)  ', 'Analytic Price  ');
    fig_num += 1;
    
endfor

plot = true;

endfunction  


%!test
%! S0 = 100;
%! r=0.0;
%! sigma=0.25;
%! T=1.1;
%! n_min = 100; % min number of time steps
%! n_max = 5000; % max number of time steps
%! delta_n = 20;
%! K_List = [85, 100, 110];
%! assert (latticeEx1(S0, r, sigma, K_List, n_min, n_max, delta_n, T))




