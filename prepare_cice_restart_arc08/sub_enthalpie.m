function qi = sub_enthalpie(fld,Ti,Tm);
% Calculates enthalipe of snow or ice (fld)
% useing formulas from CICE manual v5.1 (eq.79 - 81, 86)
% this is the amount of energy needed to melt 
% the amount of snow, for sea ice melt and 
% bring T of melted water up to 0C
% Tm - metling T, Ti - ice or snow T
c0 = 2106;    % J/(kg*Cdeg) - specific heat of fresh ice, at 0C
L0 = 3.34e5; % J/kg  latent heat of fusion of fresh ice at 0C
mu = 0.054;  % deg/ppt - liquidus, ratio betw. Tfreez and S of brine
cw = 4218;   % J/(kg*Cdeg) - specific heat of sea water

%ci = c0+L0*mu*Si/(Ti^2);
rho_snow = 330;  % kg/m3
rho_ice  = 917;  % kg/m3
%Tm       = -0.4; % ice melt T, C

switch(fld),
 case('snow'); % of fresh ice
  qi = -rho_snow*(-c0*Ti+L0); % J/m3
 case('ice');
  qi = -rho_ice*(c0*(Tm-Ti)+L0*(1-Tm/Ti)-cw*Tm);
end

% Make sure qi is <0
% usually when Tm is too low
if qi>0
  fprintf('Enthalpie > 0: Tice = %5.3f, qi = %d\n',Ti,qi);
  qi = -rho_snow*(-c0*Ti+L0); % J/m3
  fprintf('  Use snow enth, instead: %d\n',qi);
  fprintf(' ===============\n');
%  qi = -0.1;
end

% T in the layer:
% here qi is J/m3 -
% enthalpie of the whole layer
% in restart, it is divided by nlr
% and vicen = hin*aicen
% from CICE manual
a=c0;
b=(cw-c0)*Tm-qi/rho_ice-L0;
c=L0*Tm;

Tlr = (-b-sqrt(b^2-4*a*c))/(2*a);

return