clear all; close all;
graphics_toolkit("gnuplot");
pkg load financial

### k =95
S0 = 100;
K=95;
r=0.03;
sigma=0.2;
T=1.1;

n_min = 100; % min number of time steps
n_max = 5000; % max number of time steps
delta_n = 20;

[BS_C,BS_P] = blsprice(S0, K, r, T, sigma);

num_steps = linspace(n_min,n_max, ((n_max-n_min)/delta_n)+1); 
BinomPrices = zeros(1, length(num_steps));

steps = ((n_max-n_min)/delta_n)+1;


index = 1;
for n = num_steps
    delta_T = T / n;
    #BinomPrices(index) = BinEuroCall(S0, K, r, T, sigma, n); 
    BinomPrices(index) = binPriceCRR(K,S0,r,sigma,delta_T,n,'CALL',false);
    #printf("index: %d; Price: %d \n", index, BinomPrices(index));
    index = index + 1;
end;
    
figure(1)                                              
plot(num_steps, BinomPrices,'-', num_steps, ones(1,length(num_steps))*BS_C, 'b-'); 
#plot(num_steps, ones(1,length(num_steps))*BS_C, 'ro');  

### k = 100
K=100;

[BS_C,BS_P] = blsprice(S0, K, r, T, sigma);

num_steps_even = linspace(n_min,n_max, ((n_max-n_min)/delta_n)+1); 

BinomPrices_even = zeros(1, length(num_steps));
BinomPrices_odd = zeros(1, length(num_steps));
steps = ((n_max-n_min)/delta_n)+1;

index = 1;
for n_even = num_steps
    n_odd = n_even + 1;
    delta_T_even = T / n_even;
    delta_T_odd = T / n_odd;
    BinomPrices_even(index) = binPriceCRR(K,S0,r,sigma,delta_T_even, n_even,'CALL',false);
    BinomPrices_odd(index) = binPriceCRR(K,S0,r,sigma,delta_T_odd, n_odd,'CALL',false);
    index = index + 1;
end;
    
figure(2)                                              
plot(num_steps, BinomPrices_even,'-', num_steps, BinomPrices_odd, '-', num_steps, ones(1,length(num_steps))*BS_C, '-'); 

