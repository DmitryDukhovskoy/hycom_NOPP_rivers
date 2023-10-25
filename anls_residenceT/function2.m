function y = function2(x,t);
% Function 2: solution for linearly increasing GFWA
% with initial V(0) = 0 
% discharge rate
% y0 - volume of GFWA, steady-state solution (taken from HYCOM tacer experiment)
%      how much GFWA accumulated in the SPNA at time t=24 years (since 1993)
% iF0 - intercept, p - slope for regression line approximating
%       the Greenl disch anomaly - see paper for detail 
%  
iF0 = 17.44;
p   = 12.15;
y0  = 2240; % km3
y = iF0*x+p*(x*t-x^2)+(p*x^2-iF0*x)*exp(-t/x)-y0;

return



