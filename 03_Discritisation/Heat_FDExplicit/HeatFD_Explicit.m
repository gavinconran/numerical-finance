%% HeatEquationFD.m
%% Explicit finite difference equation for the heat equation on [0,1]*[0,tMax]
%% The output matrix corresponds to the solution at each grid point ((i-1)*dx, (j-1)*dt)

function heat = HeatFD_Explicit(dx, dt, tMax)
N = round(1/dx);
M = round(tMax/dt)
heat = zeros(N+1,M+1);
alpha = dt/(dx*dx)
xAxis = 0:dx:1;
for i=2:ceil((N+1)/2) %% Boundary condition at t=0
    heat(i,1) = xAxis(i)^2;
    heat(N+2-i,1) = heat(i,1);
end
% The boundary condition heat(1,j)=0=heat(N,,j) are ensured for any j

for j=1:M
    for i=2:N
        heat(i,j+1) = alpha*heat(i-1,j)+(1-2*alpha)*heat(i,j)+alpha*heat(i+1,j);
    end
end



%%***************************************
%%Plot the initial time-0 boundary condition:
%%dx=0.1;dt=0.001;heat=HeatFD_Explicit(dx, dt, 1);plot(0:dx:1,heat(:,1))

%%Plot the value function at time t=100*dt=0.1:

%dx=0.1;dt=0.001;M=round(1/dt);heat=HeatFD_Explicit(dx, dt, 1);plot(0:dx:1,heat(:,101));



