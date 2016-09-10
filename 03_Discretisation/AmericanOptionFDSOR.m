%%%%%%%%% Courtesy of Sergey Badikov (Imperial College London) %%%%%%%%%
%%%%%%%%% Date: 30 May 2014 %%%%%%%%%

%%%%%%%%% Compute the price of an American Call or Put option using a 
%%%%%%%%% Successive Over Relaxation method.
%%%%%%%%% The stock price is assumed to follow the Black-Scholes model
%%%%%%%%% without dividend.


function price = AmericanOptionFDSOR(S0,K,r,T,sigma,Smax,ds,dt,omega,tol, CallPutFlag)
%% S0: initial spot price
%% K: strike of the option
%% r: instantaneous risk-free rate
%% T: maturity of the option
%% sigma: Black-Scholes instantaneous volatility
%% Smax: maximum value for the space grid
%% ds: space step discretisation
%% dt: time step discretisation
%% omega: SOR parameter
%% tol: tolerance
%% CallPutFlag: 'True' for a Call, 'False' for a Call

%set up the grid and adjust coefficients
M = round(Smax/ds);
ds = Smax/M; 
N = round(T/dt);
dt = T/N; 

space = linspace(0,Smax,M+1)'; 
veti = 0:M;
vetj = 0:N;

% vectors for Gauss-Seidel update
oldval = zeros(M-1,1); 
newval = zeros(M-1,1) ;

%initialise the solution matix
P = zeros(M+1,N+1); 

%specify the boundary conditions
if CallPutFlag
    %Call boundary conditions
    P(:,N+1) = max(space-K,0);
    P(M+1,:) = Smax - K*exp(-r*dt*(N-vetj));
    P(1,:) = 0;
else
    %Put boundary conditions
    P(:,N+1) = max(K-space,0);
    P(1,:) = K*exp(-r*dt*(N-vetj));
    P(M+1,:) = 0;
end 

%specify the boundary conditions
payoff = max (K-space,0);
pastval = payoff; % values for the last layer
boundval = K*exp(-r*dt*(N-vetj)); % boundary values

%specify the coefficients for the Crank-Nicholson scheme 
a = 0.25*dt*( sigma^2*(veti.^2) - r*veti);
b = -dt*0.5*( sigma^2*(veti.^2) + r );
c = 0.25*dt*( sigma^2*(veti.^2) + r*veti);
M2 = diag(a(3:M) ,-1) + diag(1+b(2:M)) + diag(c(2:M-1) ,1);


x = zeros(M-1,1);
aux = zeros (M-1,1) ;

for j=N: -1 : 1
    aux(1) = a(2) * (P(1,j) + P(1,j+1));
    % set up the right-hand side and initialise
    rhs = M2*P(2:M,j+1) + aux;
    xold = P(2:M,j+1);
    error = realmax;
	% solve the sequence of linear systems by SOR method
    while tol < error
        x(1) = max ( payoff (1), xold(1) + omega/(1-b(2)) * (rhs(1) - (1-b(2))*xold(1) + c(2)*xold(2))) ;
        for k=2 : M-2
            x(k) = max ( payoff (k), xold(k) + omega/(1-b(k+1)) * (rhs(k) + a(k+1)*x(k-1) - (1-b(k+1))*xold(k) + c(k+1)*xold(k+1))) ;
        end
        x(M-1) = max( payoff (M-1), xold(M-1) + omega/(1-b(M)) * (rhs(M-1) + a(M)*x(M-2) - (1-b(M))*xold(M-1)));
        error = norm(x - xold);
        xold = x;
    end
    P(2:M,j) = x;
    x = zeros(M-1,1);
end

% return the prices, possibly by linear interpolation outside the grid
price = interp1(space, P(:,1), S0);

endfunction

%%%%%%%%%  Example
result = AmericanOptionFDSOR(50, 50, 0.1, 5/12, 0.4, 100, 2, 5/1200, 1.2, 0.001, 'True')
