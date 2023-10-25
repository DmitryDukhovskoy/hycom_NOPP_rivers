function y = function3(x,t);
% Function 3: solution for step-function increased 
% discharge rate
% initial solution V(0) = 0 
% for varying parameters
%
% y0 - volume of GFWA, steady-state solution (taken from HYCOM tacer experiment)
%      how much GFWA accumulated in the SPNA at time t=24 years (since 1993)
Phi = 166.88; % total GFWA influx into SPNA: Denm+Davis+south Greenland, km3/yr
%Phi = 175;   % GFWA = total Baffin + S greenld + fraction East Greenl
%Phi = 195; % include all Baffin Bay runoff and AO
%Phi = 209;   % total mean anomalous Greenland discharge
%y0  = 2240; % km3  
y0  = 2075; % km3 - GFWA vol from HYCOM adjusted for step-function GrFW increase
y0 = 0.8*y0;
y = Phi*x*(1-exp(-t/x))-y0;

return
