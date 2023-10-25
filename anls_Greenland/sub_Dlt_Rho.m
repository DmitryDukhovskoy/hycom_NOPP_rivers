function [Zmld,dR] = sub_Dlt_Rho(Tz,Sz,ZMz,Zz,ibb,hb,dRho0);
% Find density difference between the surface 10m layer 
% and subsurface and 
% Is bottom depth is shallower than 10m - iz10=depth of bottom
%
% MLD based on dRho0
%
% INPUT: Tz,Sz - arrays of T/S in water column with N vert. layers
%                values are given in the middle of the layers
%                or picewise continuous distribution is assumed
%                where T&S do not change within the layers
%       ZMz - array of midle layer depths (m)
%       Zz  - array of layer interface depths (m), N+1 
%       dT - T jump for defining the threshold T
%       hb - local bottom depth
%       id - grid linear index
%            for HYCOM grid (mm,nn)
%
% Depth has to be deeper than 10 m
% dT is T change wrt to Tref that defines
% MLD - density based definition
%

hb=-abs(hb);
%for diagnostics:
%[j,i]=ind2sub([mm,nn],id);


nlv=length(Tz);

% Stabilize density profiles 
% and adjust T/S profiles 
% to provide stable rho
% if needed
[Tz,Sz] = sub_stabilize_rho(Tz,Sz,ZMz,Zz);

% Assume that mixed layer is at least 10 m
iz10=max(find(ZMz>=-10));
if abs(hb)<10
  iz10=max(find(ZMz>=-abs(hb)));
end

Tref = Tz(iz10);
Sref = Sz(iz10);
%Pref = ZMz(iz10); % 1m ~ 1dBar
Rref = sw_dens0(Sref,Tref);
Rhow = sw_dens0(Sz,Tz);
dR = abs(Rref-Rhow);

%keyboard
if hb>=Zz(iz10+1),
  Zmld = hb;
  return
end

%dRho0 = 5e-3;
iz2=min(find(dR>dRho0));
if isempty(iz2);
  Zmld=hb;
  return;
end
if iz2<iz10,
  Zmld=hb;
  return;
end



iz1=iz2-1;
%keyboard
rho1 = Rhow(iz1);
rho2 = Rhow(iz2);
%keyboard
Zmld = interp1([rho1,rho2],ZMz(iz1:iz2),Rref+dRho0);

chck=0;
if chck>0
  Rho = Rhow;
  figure(10); clf;
  axes('position',[0.08 0.08 0.6 0.82]);
  plot(Rhow,ZMz,'.-');
  hold;
  plot([Rref+dRho0 Rref+dRho0],[min(ZMz) max(ZMz)],'r--');
  plot([Rref Rref],[min(ZMz) max(ZMz)],'m:');
  set(gca,'ylim',[hb 0]);
  title('MLD(g), RhoRef(magenta), RhoMLD(r)');
  btx='sub_mld_drho.m';
  bottom_text(btx,'pwd',1);
end

dR(ibb+1:end)=nan;


return