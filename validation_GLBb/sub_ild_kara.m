function Zild = sub_ild_kara(Tz,Sz,ZMz,Zz,dT,hb,id,mm,nn);
% Follow methodology from Kara et al. (2000)
% Depth has to be deeper than 10 m
% dT is T change wrt to Tref that defines
% Isothermal Layer Depth (ILD)
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
% Note that Kara's algorithm for IDL does not
% work well for polar regions due to small T varia
% tions and T can have several extremum (e.g. Atl water +, then -)
% within a relative dltT range
% need to keep dT small (<0.5)
%
% D. Dukhovskoy, COAPS FSU, June 2017
%

%for diagnostics:
[j,i]=ind2sub([mm,nn],id);


% Stabilize density profiles 
% and adjust T/S profiles 
% to provide stable rho
% if needed
[Tz,Sz] = sub_stabilize_rho(Tz,Sz,ZMz,Zz);

% Step 1: Define Tref
% Update Tref = T (z) as move
% downward
iz10=max(find(ZMz>-11));
Tref = Tz(iz10);

dlt = 0.1*dT;   % reference difference

nlv=length(Tz);

% If too shallow - skip
if hb>=Zz(iz10+1),
  Zild = hb;
  return
end


% Update Tref
% by looking for a uniform T region:
for kk=iz10+1:nlv
  dt0 = abs(Tz(kk-1)-Tz(kk));
  sdt = Tz(kk-1)-Tz(kk);
  
  if dt0<=dlt % homogeneous Layer
    Tref=Tz(kk-1); % update ref T
  else
    break
  end
  
end

% Find dt/dz below the homogeneous layer >0/<0
%tmm = nanmean(Tz(kk:end));
%sdt = (Tz(kk-1)-tmm);
%sdt = -nanmean(diff(Tz(kk-1:end)));

% Define T of the isotherm. mixed layer:
% This approach described in Kara et al 2000
% does not work for cases of small dT jumps along T profile
% as T(kk) can be > T(kk-1) but below kk, all T may decrease
% and Tb will be wrongly defined (as it is Tb=Tref+dt , whereas
% should be Tref-dt).
% Instead, look for abs. difference dT from Tref
%if sdt<0 % T increases with depth
%  Tb = Tref-dT; % 
%else % T decreases with depth
%  Tb = Tref+dT; % 
%end

% Find depth of the Isotherm. Layer
% Locate layer where abs(dltT) exeeds dT
iz2=[];
for iz=kk:nlv
  if isnan(Tz(iz)); continue; end
  dt0=Tz(iz)-Tref;
  if abs(dt0)>dT
    iz2=iz;
    if dt0<0,
      Tb = Tref-dT;
    else
      Tb = Tref+dT;
    end
    break;
  end
end

% Use Tref at 10m
% if no homogeneous layer found below 10m
if isempty(iz2) 
  Tref = Tz(iz10);
  for iz=iz10+1:nlv
    if isnan(Tz(iz)); continue; end
    dt0 = Tz(iz)-Tref;
    
    if abs(dt0)>dT
      iz2=iz;
      if dt0<0,
	Tb = Tref-dT;
      else
	Tb = Tref+dT;
      end
      break;
    end
  end

% if no dt0>dT found and hit the bottom, ILD = total depth
%keyboard
  if abs(dt0)<dT
%    Zild = ZMz(nlv);
    fprintf('!!! i=%i,j=%i No isothermal L. found: ILD to bottom %5.1f\n',i,j,hb);
%    if i>600 | j>1000
%    keyboard; % check profile
%    end;
    Zild=hb;
    return;
  end
end
iz1=iz2-1;
%keyboard
Zild = interp1(Tz(iz1:iz2),ZMz(iz1:iz2),Tb);

chck=0;
if chck>0
  figure(10); clf;
  axes('position',[0.08 0.08 0.6 0.82]);
  plot(Tz,ZMz,'.-');
  hold;
  plot([Tb Tb],[min(ZMz) max(ZMz)],'r--');
  plot([min(Tz) max(Tz)],[Zild Zild],'g--');
  plot([Tref Tref],[min(ZMz) max(ZMz)],'m:');
  set(gca,'ylim',[hb 0]);
  title('ILD(g), TRef(magenta), T-MLD Base(r)');
  btx='sub_ild_kara.m';
  bottom_text(btx,'pwd',1);
end



return