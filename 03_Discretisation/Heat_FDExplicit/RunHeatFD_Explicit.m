%% Temp_Jack.m
clear all; close all;
graphics_toolkit("gnuplot");
dx=0.1;dt=0.004;tMax=1;
%sol=HeatFD_Explicit(dx, dt, tMax);
%subplot(2,2,1);
%plot(0:dx:1,sol(:, 1))
%subplot(2,2,2);
%plot(0:dx:1,sol(:,11))
%subplot(2,2,3);
%plot(0:dx:1,sol(:,51))
%subplot(2,2,4);
%plot(0:dx:1,sol(:,101))

sol=HeatFD_Explicit(dx, dt, tMax);
subplot(2,1,1);
plot(0:dx:1,sol(:, 1))
subplot(2,1,2);
plot(0:dx:1,sol(:,size(sol,2)))
