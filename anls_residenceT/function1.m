function y = function1(x,t);
% Function 1: solution for step-function increased 
% discharge rate
% initial solution V(0) = 0 
%
% y0 - volume of GFWA, steady-state solution (taken from HYCOM tacer experiment)
%      how much GFWA accumulated in the SPNA at time t=24 years (since 1993)
Phi = 166.88; % total GFWA influx into SPNA: Denm+Davis+south Greenland, km3/yr
%y0  = 2240; % km3  
y0  = 2075; % km3 - GFWA vol from HYCOM adjusted for step-function GrFW increase
y = Phi*x*(1-exp(-t/x))-y0;

return
