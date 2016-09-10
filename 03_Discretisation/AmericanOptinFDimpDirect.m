%%%%%%%%% Courtesy of Sergey Badikov (Imperial College London) %%%%%%%%%
%%%%%%%%% Date: 30 May 2014 %%%%%%%%%

%%%%%%%%% Compute the price of an American Call or Put option using an implicit
%%%%%%%%% finite-difference scheme.
%%%%%%%%% The stock price is assumed to follow the Black-Scholes model
%%%%%%%%% without dividend.

function price = AmericanOptinFDimpDirect(S0, K, r, T, sigma, Smax, ds, dt, CallPutFlag)
%% S0: initial spot price
%% K: strike of the option
%% r: instantaneous risk-free rate
%% T: maturity of the option
%% sigma: Black-Scholes instantaneous volatility
%% Smax: maximum value for the space grid
%% ds: space step discretisation
%% dt: time step discretisation
%% CallPutFlag: 'True' for a Call, 'False' for a Call


%set up the grid and adjust coefficients
M = round(Smax/ds); %% space discretisation grid
M
ds = Smax/M; 
N = round(T/dt); %% time discretisation grid
N

dt = T/N; 

space = linspace(0,Smax,M+1)'; 

veti = 0:M;
vetj = 0:N;

%initialise the solution matrix
P = zeros(M+1,N+1); 

%initialise the boundary conditions
if CallPutFlag
    %call boundary conditions
    P(:,N+1) = max(space-K,0);
    P(M+1,:) = Smax - K*exp(-r*dt*(N-vetj));
    P(1,:) = 0;
else
    %put boundary conditions
    P(:,N+1) = max(K-space,0);
    P(1,:) = K*exp(-r*dt*(N-vetj));
    P(M+1,:) = 0;
end 

%define the coefficients for the iteration matrix of the implicit scheme 
a = 0.5*dt*(r*veti - sigma^2*(veti.^2)); 
b = 1 + dt*(sigma^2*veti.^2 + r);
c = -0.5*dt*(sigma^2*veti + r).*veti;

%define the iteration matrix 
A = diag(a(3:M),-1) + diag(b(2:M)) + diag(c(2:M-1),1);

aux = zeros(M-1,1);
[L, U] = lu(A); 

%%%%%%%%%%%% run the the implicit scheme backwards from the maturity, and
%%%%%%%%%%%% compare at each time step with the exercice price.

for j=N:-1:1
   aux(1) =  a(2) * P(1,j);
   aux(M-1) = c(M) * P(M+1,j);
   P(2:M,j) = U \ (L \ (P(2:M, j+1) - aux)); 
   for i=2:M
       P(i,j) = max(P(i,j),P(i,N+1));
   end

end

%%%%%%%%%%%% returns the price of the American option
price = interp1(space, P(:,1), S0);
size(P)

endfunction

%%%%%%%%%%%% Example %%%%%%%%%%%%
result = AmericanOptinFDimpDirect(50, 50, 0.1, 5/12, 0.4, 100, 2, 5/1200,  'True');
result

