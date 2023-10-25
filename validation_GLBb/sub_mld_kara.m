function Zmld = sub_mld_kara(Tz,Sz,ZMz,Zz,dT,hb,id,mm,nn);
% Follow methodology from Kara et al. (2000)
% Note that in Kara's approach, d(RHo) is T based
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
% Isothermal Layer Deith (ILD) - see sub_ild_kara.m
% MLD - density based definition
%  rho(T+dt,S,P)-rho(T,S,P)
%  See Kara et al., 2000 - state that 
%  as T is potential in most obs., use P=0
%  incompressibility assumption for rho agrees
%  with most models
%

%for diagnostics:
[j,i]=ind2sub([mm,nn],id);

% Step 1: Define Tref
% Update Tref = T (z) as move
% downward
nlv=length(Tz);

% Stabilize density profiles 
% and adjust T/S profiles 
% to provide stable rho
% if needed
[Tz,Sz] = sub_stabilize_rho(Tz,Sz,ZMz,Zz);

iz10=max(find(ZMz>=-10));
Tref = Tz(iz10);
Sref = Sz(iz10);
%Pref = ZMz(iz10); % 1m ~ 1dBar
Rref = sw_dens0(Sref,Tref);
dR = abs(Rref-sw_dens0(Sref,Tref-dT));

if hb>=Zz(iz10+1),
  Zmld = hb;
  return
end

dlt = 0.1*dT;   % reference difference
dRho0 = abs(sw_dens0(Sref,Tref)-sw_dens0(Sref,Tref-dlt));

% Update Tref
% by looking for a uniform T region:
for kk=iz10+1:nlv
  t1=Tz(kk-1);
  t2=Tz(kk);
  s1=Sz(kk-1);
  s2=Sz(kk);
  drho = abs(sw_dens0(s2,t2)-sw_dens0(s1,t1));
  if drho<dRho0, % homogen. layer
    Tref=Tz(kk-1); % update ref T
    Sref=Sz(kk-1); % update ref S
    Rref = sw_dens0(Sref,Tref);% ref dens.
    dR = abs(Rref-sw_dens0(Sref,Tref-dT));
  else
    break
  end
    
end


% Find depth of the Isotherm. Layer
% Locate layer where abs(dltT) exceeds dR
iz2=[];
for iz=kk:nlv
  if isnan(Tz(iz)); continue; end
  t2=Tz(iz);
  s2=Sz(iz);
  dr0 = abs(sw_dens0(s2,t2)-Rref);
  
  if dr0>dR
    iz2=iz;
    Rb = sw_dens0(Sref,Tref-dT);
    break;
  end
end

% Use Tref at 10m
% if no homogeneous layer found below 10m
if isempty(iz2) 
  Rref = sw_dens0(Sz(iz10),Tz(iz10));
  for iz=iz10+1:nlv
    if isnan(Tz(iz)); 
      dr0=dR-0.0001;
      iz2=[];
      break; % nan values = hit bottom 
    end
    t2=Tz(iz);
    s2=Sz(iz);
    dr0 = abs(sw_dens0(s2,t2)-Rref);
    
    if dr0>dR
      iz2=iz;
      Rb = sw_dens0(Sref,Tref-dT);
      break;
    end
    
  end

%  if ~exist('dr0','var'); keyboard; end;
% if no dr0>dR found and hit the bottom, ILD = total depth
  if (dr0<dR | isempty(iz2)) 
%    Zmld = ZMz(nlv);
    if (abs(hb)>1000)
    fprintf('*** i=%i, j=%i: No picnocline found: MLD = bottom %5.1f\n',i,j,hb);
    end
%    keyboard; % check profile
    Zmld=hb;
    return;
  end
end
iz1=iz2-1;

rho1 = sw_dens0(Sz(iz1),Tz(iz1));
rho2 = sw_dens0(Sz(iz2),Tz(iz2));
%keyboard
Zmld = interp1([rho1,rho2],ZMz(iz1:iz2),Rb);

chck=0;
if chck>0
  Rho = sw_dens0(Sz,Tz);
  figure(10); clf;
  axes('position',[0.08 0.08 0.6 0.82]);
  plot(Rho,ZMz,'.-');
  hold;
  plot([Rb Rb],[min(ZMz) max(ZMz)],'r--');
  plot([min(Rho) max(Rho)],[Zmld Zmld],'g--');
  plot([Rref Rref],[min(ZMz) max(ZMz)],'m:');
  set(gca,'ylim',[hb 0]);
  title('MLD(g), RhoRef(magenta), RhoMLD Base(r)');
  btx='sub_mld_kara.m';
  bottom_text(btx,'pwd',1);
end



return