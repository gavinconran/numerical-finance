% 1Lattice. Exercise 5 
% Exercise 5. We plot here twenty Black-Scholes paths, where the time interval [0, 1] 
% has been split respectively into ten, one hundred and one thousand subintervals.

clear all; close all;
graphics_toolkit("gnuplot");
pkg load financial

function plot = MonteCarloBSPath(S0, r, sigma, Tsubs_List, T, sims)

% Inputs : S0 - stock price
%        : r - risk free interest rate
%        : sigma - volatility
%        : Tsubs_List - List of number of Time sub intervals
%        : T - time duration
%
% Output: plot - true or false


max_steps = Tsubs_List(end);

fig_num = 1; % update figure number for plotting
for steps = Tsubs_List
    % Create a vector of number of steps
    time_steps = linspace(0, 1, steps+1); 
    % Initialise the vector of Option Prices
    Price_t = zeros(sims, length(time_steps));
    Price_t(1,:) = S0;
    t = T / steps;;
    
    % Plot    
    figure(fig_num)     
    for i=1:sims
        % For each time step compute the Stock Price St 
        for j = 2:length(time_steps)            
            Price_t(i,j) = Price_t(i,j-1) + sqrt(t)*normrnd(r,sigma); 
        endfor
    endfor 
                                             
    plot(time_steps, Price_t,'-');
    xlabel('Time (T)');
    ylabel('Stoke Price');
    str = sprintf('Black Scholes path for Stock Price. Sub Interval = %d', steps);
    title(str); 
    ylim([-3, 3]);
    fig_num += 1;
endfor

plot = true;

endfunction  

%!test
%! S0 = 0;
%! r=0.0;
%! sigma=0.10;
%! T=1.0;
%! Tsubs_List = [10, 100, 1000];
%! sims = 20;
%! assert (MonteCarloBSPath(S0, r, sigma, Tsubs_List, T, sims))




