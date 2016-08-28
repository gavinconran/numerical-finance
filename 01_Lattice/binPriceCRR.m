%%% obtained from http://www.goddardconsulting.ca/matlab-binomial-crr.html

function oPrice = binPriceCRR(X,S0,r,sig,dt,steps,oType,earlyExercise)
% Function to calculate the price of a vanilla European or American
% Put or Call option using a Cox Ross Rubinstein binomial tree.
%
% Inputs: X - strike
%       : S0 - stock price
%       : r - risk free interest rate
%       : sig - volatility
%       : dt - size of time steps
%       : steps - number of time steps to calculate
%       : oType - must be 'PUT' or 'CALL'.
%       : earlyExercise - true for American, false for European.
%
% Output: oPrice - the option price
%
% Notes: This code focuses on details of the implementation of the 
%                               Cox Ross Rubinstein (CRR)
%        algorithm.
%        It does not contain any programatic essentials such as error
%        checking.
%        It does not allow for optional/default input arguments.
%        It is not optimized for memory efficiency or speed.

% Author: Phil Goddard (phil@goddardconsulting.ca)
% Date  : Q4, 2007

% Calculate the Cox Ross Rubinstein model parameters
a = exp(r*dt);          % discount rate
u = exp(sig*sqrt(dt));  % up
d = 1/u;                % down
p = (a-d)/(u-d);        % martingale probability

% Loop over each node and calculate the Cox Ross Rubinstein underlying price tree
priceTree = nan(steps+1,steps+1);
priceTree(1,1) = S0;
for idx = 2:steps+1
    priceTree(1:idx-1,idx) = priceTree(1:idx-1,idx-1)*u;
    priceTree(idx,idx) = priceTree(idx-1,idx-1)*d;
end

% Calculate the value at expiry
valueTree = nan(size(priceTree));
switch oType
    case 'PUT'
        valueTree(:,end) = max(X-priceTree(:,end),0);
    case 'CALL'
        valueTree(:,end) = max(priceTree(:,end)-X,0);
end

% Loop backwards to get values at the earlier times
steps = size(priceTree,2)-1;
for idx = steps:-1:1
    valueTree(1:idx,idx) = ...
        exp(-r*dt)*(p*valueTree(1:idx,idx+1) ...
        + (1-p)*valueTree(2:idx+1,idx+1));
    if earlyExercise
        switch oType
            case 'PUT'
                valueTree(1:idx,idx) = ...
                    max(X-priceTree(1:idx,idx),valueTree(1:idx,idx));
            case 'CALL'
                valueTree(1:idx,idx) = ...
                    max(priceTree(1:idx,idx)-X,valueTree(1:idx,idx));
        end
    end
end

% Output the option price
oPrice = valueTree(1);   
