function Zmld = sub_mld_DD(Tz,Sz,ZMz,Zz,hb);
% Find density gradient
% Fit polynomial and find stDev
% Upper M.L - is where dRho/dz < 
% 
% D. Dukhovskoy, COAPS FSU, June 2017
%
% downward
iz10=max(find(ZMz>=-10))+1;
Tref = Tz(iz10);
Sref = Sz(iz10);

nlv=length(Tz);

% Stabilize density profiles 
% and adjust T/S profiles 
% to provide stable rho
% if needed
%keyboard
Tz0=Tz;
Sz0=Sz;
[Tz,Sz] = sub_stabilize_rho(Tz,Sz,ZMz,Zz);

% Calculate dRho/dz: fwd Euler difference,
% backward Euler difference
% interpolate into pnt ik (averge two)
% 
Rho = sw_dens0(Sz,Tz);
dRdz = Tz*nan;
dRdz(1:iz10)=0;
for ik=iz10+1:nlv-1
%  fprintf('ik=%i\n',ik);
% Fwd Euler:
  dzf = abs(ZMz(ik+1)-ZMz(ik));
  dRf = (Rho(ik+1)-Rho(ik))/dzf;
  if isnan(dRf)| dzf<0.01; 
    if ~isnan(Tz(ik+1))
      dRdz(ik) = dRdz(ik-1); % near the bottom
    end
    break; 
  end;
  zf  = 0.5*abs(ZMz(ik+1)+ZMz(ik)); % depth of Fwd derivative
  dzb = abs(ZMz(ik)-ZMz(ik-1));
  dRb = (Rho(ik)-Rho(ik-1))/dzb;
  zb  = 0.5*abs(ZMz(ik-1)+ZMz(ik));
  dRz = dRb+(dRf-dRb)/(zf-zb)*(abs(ZMz(ik))-zb); % interpolate into depth ZMz
  if sign(dRz)~=sign(0.5*(dRb+dRf))
    fprintf(' dRdz <0 , %6.3d, Making it 1e-6  Layer %i!!!\n',dRz,ik);
    dRz=1e-9;
%    keyboard
%    error('check dRdz calculation - sign derivative wrong!!!');
  end
  
  dRdz(ik) = dRz;
end

%Nfr = sqrt(9.81/1026*dRdz); % Brunt-Viasala freq.

% overall Mean gradient & stDev:
dRmn = nanmean(dRdz);
stR  = nanstd(dRdz);
stR10= nanstd(dRdz(1:iz10)); 
r25  = prctile(dRdz,25);
r50  = prctile(dRdz,50);
r75  = prctile(dRdz,75);
iqr  = abs(r75-r25);

% Fix several cases of 0 grad Rho - errors
% in T/S fields resulted in inverse Rho
% that after stabilization have completely mixed 
% water column
if stR<1e-12; 
  stR = 1e-12;
end
if stR10<1e-12; 
  stR10 = 1e-12;
end
%keyboard
% Assume that strong pycnocline > 2 stDev dRho/dz:
% If not such value than the water column is close to stably stratified
% Find depth of dRdz0
dRdz0 = dRdz(iz10)+2*stR; % seems to work the best
%dRdz0 = dRdz(iz10)+1.5*iqr;
%
% Alternatively: Look for values exceeding 1.5 interquartile range
% from the 75th percentile
%dRdz0 = r75+0.5*iqr; % too extreme, especially if distr is too skewed -> high iqr


% Find depth of the reference dRho/dz
dRz = 0;
iz2 = [];
for iz=iz10+1:nlv
  if isnan(dRdz(iz)); continue; end
  dRz = dRdz(iz);
  if dRz>=dRdz0
    iz2=iz;
    break;
  end
end

% if no dRz>dRdz0 found 
% MLD = depth where dRz > 2*std Rho(upper 10 m)
if (dRz<dRdz0 | isempty(iz2)) 
   if abs(hb)>1000
%    Zmld = ZMz(nlv);
    fprintf('xxx DD:: No picnocl found, max grad(Rho)= %8.6f MLD=%6.2f\n',...
	  max(dRdz),hb);  
  end
  
  dRdz10 = dRdz(iz10)+2*stR10;
  iz2 = min(find(dRdz>=dRdz10));
  if isempty(iz2)
    Zmld=hb;
    return;
  end
  
%  keyboard; % check profile
  dRdz0 = dRdz10;
  
end


iz1=iz2-1;
drho1 = dRdz(iz1); 
drho2 = dRdz(iz2);
Zmld = interp1([drho1,drho2],ZMz(iz1:iz2),dRdz0);

chck=0;
if chck>0
%  Rho = sw_dens0(Sz,Tz);
  figure(10); clf;
  axes('position',[0.08 0.08 0.6 0.82]);
  plot(Rho,ZMz,'.-');
  hold;
%  plot([Rb Rb],[min(ZMz) max(ZMz)],'r--');
  plot([min(Rho) max(Rho)],[Zmld Zmld],'g--');
%  plot([Rref Rref],[min(ZMz) max(ZMz)],'m:');
  set(gca,'ylim',[hb 0]);
  title('MLD(g), RhoRef(magenta), RhoMLD Base(r)');
  btx='sub_mld_DD.m';
  bottom_text(btx,'pwd',1);
  keyboard
end



return