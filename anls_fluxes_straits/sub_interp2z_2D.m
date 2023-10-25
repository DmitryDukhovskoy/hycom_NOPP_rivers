function Fi = sub_sct_interp2z(SCT,F, ZZf,ZMf,ZZh,ZMh);
%
% Interpolate normal to a vertical section 
% (2D array) from HYCOM, hybr grid
% SCT - structured aray with info about the section
%
% onto z-fixed depths (ZZf)
% ZZh - HYCOM 2D array of layer interface depths
% ZMh - HYCOM 2D array of rho-point depths
% ZZf - 1D array of depths fixed depths 
% ZMf - 1D array of mid-layer depths for fixed layers
%
% Interpolated field F is in the middle point
% of the segments - ready for calculating 
% fluxes, thus dim of F is -1 compared to 
% the number of segment points
%
%fprintf('Interpolating U,T,S hybrid grid -> z\n');

hg = 1e20;
[nlr,nn] = size(F);
nZ = length(ZZf);
mZ = length(ZZf)-1;
Zmin=min(ZZf);
%keyboard

Hs = SCT.SGM_Hb;
np = length(Hs);

Xf=[1:nn];
%[XXi,ZMi]=meshgrid(Xf,ZMf);

% Add extra bottom layer for interpolation
% below the min depth on fixed Z-array
ZMh(nlr+1,:)=-12000;
F(nlr+1,:)=F(nlr,:);

% Add extra surface layer for interpolation
ZMh=[zeros(1,np);ZMh];
F=[F(1,:);F];

In=find(isnan(F));
if ~isempty(In)
  F=sub_fill_bottom_nans(F);
%  fprintf(' sub_inter2z_2D: *** ERR *** found nans in the field\n');
%  keyboard
%  error(' STOPPING INTERPOLATION');
end

[nlr,nn] = size(F);
[Xh,dmm]=meshgrid(Xf,[1:nlr]);



% Interpolation:
for kk=1:np
  zm=ZMh(:,kk);
% For interpolation
% depths cannot be the same
  dz=abs(diff(zm));
  Iz=find(dz==0);
  if ~isempty(Iz)
    nz=length(Iz);
    for kz=1:nz
      iz0=Iz(kz)+1;
      zm(iz0)=zm(iz0-1)-0.001;
    end
  end
  
  f0=F(:,kk);
  fi=interp1(zm,f0,ZMf,'pchip');
  Fi(:,kk)=fi;
end


return