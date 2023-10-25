function FWF = sub_calc_sflx(S,U,V,dH,SCT,HH,Sref,DX,DY,z0);
% integrate to z0 and to bottom
% S or any other tracer
% but change Sref = 0 if not salinity
%
z0 = abs(z0);
[nlr,mm,nn] = size(U);
II = SCT.II;
JJ = SCT.JJ;
Hs = SCT.Hbottom;
np = length(II);
Vflx = zeros(1,np); % volume flux for checking
FWF.VolFlx  = zeros(1,np);
FWF.VolFlxZ = zeros(1,np);
FWF.FWFlx   = zeros(1,np);
FWF.FWFlxZ  = zeros(1,np);

for ig=1:np-1
  i1=II(ig);
  i2=II(ig+1);
  j1=JJ(ig);
  j2=JJ(ig+1);

  if HH(j1,i1)>=0 & HH(j2,i2)>=0, 
    continue; 
  end;

  if Hs(ig)>0, continue; end;

  xsct = 1;
  if i2==i1 % y section
    xsct = 0;
  end
% Reorder i1/j1 and i2/j2 so that i2>i1 or j1>j2
  if i2<i1, 
    imm=i2;
    i2=i1;
    i1=imm;
  end
  if j2<j1
    jmm=j2;
    j2=j1;
    j1=jmm;
  end

  nrm = SCT.Norm(ig,:);
  nx  = nrm(1);
  ny  = nrm(2);
  dZ = 0;

% Collocate U, V with tracers
% Calculate fluxes
% then interpolated into middle of segment (i,i+1)
%
%                 T,dH,S,Tr 
%      - u(i,j2)   * (i,j2)    - u(i+1,j2)
%
%            <---- x  Flux at this location, interpolate from T-pnt (j+1) and (j) 
%
%      - u(i,j1)   * (i,j1)    - u(i+1,j)
%
  if xsct==0; % y-section, i1=i2
% Point (i,j1):    
    u10   = squeeze(U(:,j1,i1))*nx;
    dh10  = 0.5*squeeze(dH(:,j1,i1-1)+dH(:,j1,i1));
    u1p1  = squeeze(U(:,j1,i1+1))*nx;
    dh1p1 = 0.5*squeeze(dH(:,j1,i1)+dH(:,j1,i1+1));
    u11   = (u10.*dh10+u1p1.*dh1p1)./(dh10+dh1p1); 
    u11(dh10==0 & dh1p1==0) = nan;
% Point (i,j2):    
    u20   = squeeze(U(:,j2,i2))*nx;    
    dh20  = 0.5*squeeze(dH(:,j2,i2-1)+dH(:,j2,i2));
    u2p1  = squeeze(U(:,j2,i2+1))*nx;    
    dh2p1 = 0.5*squeeze(dH(:,j2,i2)+dH(:,j2,i2+1));
    u22   = (u20.*dh20+u2p1.*dh2p1)./(dh20+dh2p1);    
    u22(dh20==0 & dh2p1==0) = nan;
% Depth-integrate and
% Interpolate flux into mid-pnts of section (j:j+1):
    dy12 = DY(j1,i1);
    dh11 = squeeze(dH(:,j1,i1));
    dh22 = squeeze(dH(:,j2,i2));
    cdh1 = cumsum(dh11);
    iz0  = max(find(cdh1<=z0));
% Vol flux: m3/s    
    Vf   = nansum(0.5*(u11.*dh11+u22.*dh22)*dy12);
    Vfz  = nansum(0.5*(u11(1:iz0).*dh11(1:iz0)+...
		       u22(1:iz0).*dh22(1:iz0))*dy12);
% FW flux, m3/s    
    s11  = squeeze(S(:,j1,i1));
    s22  = squeeze(S(:,j2,i2));
    if Sref>0 % S
      s11(s11>Sref) = nan;
      s22(s22>Sref) = nan;
      Ff11 = (Sref-s11)/Sref.*u11.*dh11; %FWflux
      Ff22 = (Sref-s22)/Sref.*u22.*dh22;
    else  % not S
      Ff11 = s11.*u11.*dh11; %Flx of C (kg/m3) = kg/s*m
      Ff22 = s22.*u22.*dh22; %Flx of C (kg/m3) = kg/s*m
    end      
    Ff   = nansum(0.5*(Ff11+Ff22)).*dy12; % m3/s for S or kg/s
    Ffz  = nansum(0.5*(Ff11(1:iz0)+Ff22(1:iz0))).*dy12;
    
  else  % x-section, j1=j2  
%
%     *(j+1,i-1)        *(j+1,i)      *(j+1,i+1)
%
%
%     v(j+1,i-1)      v(j+1,i)     v(j+1,i+1)     
%     |                 |       ^     |
%                               | 
%     *(j,i-1)          *(j,i)  x     *(j,i+1)
%
%     |                 |             |
%     v(j+1,i-1)      v(j+1,i)     v(j+1,i+1)     
% Point (i1,j1):    
    v10   = squeeze(V(:,j1,i1))*ny;
    dh10  = 0.5*squeeze(dH(:,j1-1,i1)+dH(:,j1,i1));
    v1p1  = squeeze(V(:,j1+1,i1))*ny;
    dh1p1 = 0.5*squeeze(dH(:,j1,i1)+dH(:,j1+1,i1));
    v11   = (v10.*dh10+v1p1.*dh1p1)./(dh10+dh1p1); 
    v11(dh10==0 & dh1p1==0) = nan;

% Point (i+1,j)    
    v20   = squeeze(V(:,j1,i1+1))*ny;
    dh20  = 0.5*squeeze(dH(:,j1-1,i1+1)+dH(:,j1,i1+1));
    v2p1  = squeeze(V(:,j1+1,i1+1))*ny;
    dh2p1 = 0.5*squeeze(dH(:,j1,i1+1)+dH(:,j1+1,i1+1));
    v22   = (v20.*dh20+v2p1.*dh2p1)./(dh20+dh2p1); 
    v22(dh10==0 & dh1p1==0) = nan;

% Depth-integrate and
% Interpolate flux into mid-pnts of section (j:j+1):
    dx12 = DX(j1,i1);
    dh11 = squeeze(dH(:,j1,i1));
    dh22 = squeeze(dH(:,j2,i2));
    cdh1 = cumsum(dh11);
    iz0  = max(find(cdh1<=z0));
%
% Vol Flux, m3/s
    Vf   = nansum(0.5*(v11.*dh11+v22.*dh22)*dx12);
    Vfz  = nansum(0.5*(v11(1:iz0).*dh11(1:iz0)+...
		       v22(1:iz0).*dh22(1:iz0))*dx12);
  
% FW Flux, m3/s    
    s11  = squeeze(S(:,j1,i1));
    s22  = squeeze(S(:,j2,i2));
    if Sref>0 % S
      s11(s11>Sref) = nan;
      s22(s22>Sref) = nan;
      Ff11 = (Sref-s11)/Sref.*v11.*dh11; %FWflux
      Ff22 = (Sref-s22)/Sref.*v22.*dh22;
    else
      Ff11 = s11.*v11.*dh11; %Tracer flx kg/s*m
      Ff22 = s22.*v22.*dh22;
    end      
    Ff   = nansum(0.5*(Ff11+Ff22)).*dx12;
    Ffz  = nansum(0.5*(Ff11(1:iz0)+Ff22(1:iz0))).*dx12; % m3/s or kg/s
  end  % if
  
  FWF.VolFlx(ig)  = Vf;
  FWF.VolFlxZ(ig) = Vfz;
  FWF.FWFlx(ig)   = Ff;
  FWF.FWFlxZ(ig)  = Ffz;
end % for ig - segments

    
return