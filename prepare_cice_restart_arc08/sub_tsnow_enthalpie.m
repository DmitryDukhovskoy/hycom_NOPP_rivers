function [Tsn,Tmax] = sub_tsnow_enthalpie(esnon,vsnon,aicen);
% Follow CICE to get Tsnow from enthalpie
%
%
cp_ice = 2106;    % J/(kg*Cdeg) - specific heat of fresh ice, at 0C
L0 = 3.34e5; % J/kg  latent heat of fusion of fresh ice at 0C
mu = 0.054;  % deg/ppt - liquidus, ratio betw. Tfreez and S of brine
cw = 4218;   % J/(kg*Cdeg) - specific heat of sea water
%ci = c0+L0*mu*Si/(Ti^2);
rho_snow = 330;  % kg/m3
rho_ice  = 917;  % kg/m3
puny = 1e-11; 

nslyr = 1; % # of snow layers
hs_min = 1e-4;
hsn = vsnon/aicen; % snow thickness
hslyr = hsn/nslyr; 

if hslyr<hs_min; 
  qsn = -L0*rho_snow;
else
  qsn = esnon./vsnon;
end

Tmax = -qsn*puny*nslyr/(rho_snow*cp_ice*vsnon); % possible max T for qsn
Tsn = (L0 + qsn/rho_snow)/cp_ice;
%keyboard

return