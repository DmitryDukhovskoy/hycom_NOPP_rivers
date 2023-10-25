function FSGM = sub_vhFlx_sgm(GC,U,V,S,T,dH,DX,DY,HH,VCT,Tref);
% Calculate vol and heat fluxes across segments
% of the contour
Cp = 4200; % J/kg K
Vct = max(max(abs(VCT)));
fprintf('=== Adjusting biases in fluxes across segments,max|Vct|=%12.10f ...\n',Vct);

[nlr,mm,nn] = size(U);

II = GC.cntr_Iindx;
JJ = GC.cntr_Jindx;
Hs = GC.Hbottom;

np = length(II);
Hflx = zeros(nlr,np)*nan;
ZZ = zeros(nlr,np)*nan;
Asgm = zeros(nlr,np);
Vflx = zeros(1,np); % volume flux for checking
for ig=1:np-1
  i1=II(ig);
  i2=II(ig+1);
  j1=JJ(ig);
  j2=JJ(ig+1);

  clear vflx u12 u34 v12 v34 

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

  nrm = GC.Norm_in(ig,:);
  nx  = nrm(1);
  ny  = nrm(2);
  dZ = 0;
  clear hf* vf*
  Vct = VCT(:,ig);

% Fluxes calculated:
%            v(i,j+1)        v(i+1,j+1)
%      |--------|-------|-------|---------|
%      |                |                 |
%      |                |                 |
%      |     V1         |       V2        |  Flx = Cp*rho*(T-Tref)V1*dH1*dX1/2+
%      -        *===============*         |        Cp*rho*(T-Tref)V2*dH2*dX2/2+
%u(i,j)|     T,S,dH     |                 |   where V1 and V2 are flux-averaged
%      |                |                 |   v(i,j+1)&v(i,j) and v(i+1,j+1)&v(i+1,j)
%      |                |                 |
%      |--------|-------|-------|---------|
%      |      v(i,j)         v(i+1,j)
%      |   

  if xsct==0 % i2=i1
    if nx==0, error('Check norm x comp., it is 0!'); end;
% U Fluxes in grid cell (j1,i1):	
    u1 = squeeze(U(:,j1,i1))*nx;
    dh1= 0.5*squeeze(dH(:,j1,i1-1)+dH(:,j1,i1));
    u1(dh1<1e-3)=nan;
    u1  = sub_chck_uv(u1,dh1,0);

    u2 = squeeze(U(:,j1,i1+1))*nx;
    dh2= 0.5*squeeze(dH(:,j1,i1)+dH(:,j1,i1+1));
    u2(dh2<1e-3) = nan;
    u2  = sub_chck_uv(u2,dh2,0);	

% Values in the center of the grid cell	
    u12 = (u1.*dh1+u2.*dh2)./(dh1+dh2);
    u12 = u12+Vct; % correction for biased U to get ~=0 overall transp
    dh12= squeeze(dH(:,j1,i1));
    t12 = squeeze(T(:,j1,i1));
    s12 = squeeze(S(:,j1,i1));
    rho12= sw_dens0(s12,t12);
    dy12 = DY(j1,i1);
    hf12 = Cp*rho12.*(t12-Tref).*u12.*dh12*dy12; % W=J/s
    vf12 = u12.*dh12.*dy12;                    % m3/s - volume flux

% U Fluxes in grid cell (j2,i2):	
    u3  = squeeze(U(:,j2,i2))*nx;
    dh3 = 0.5*squeeze(dH(:,j2,i2-1)+dH(:,j2,i2));
    u3(dh3<1e-3) = nan;
    u3  = sub_chck_uv(u3,dh3,0);

    u4  = squeeze(U(:,j2,i2+1))*nx;
    dh4 = 0.5*squeeze(dH(:,j2,i2)+dH(:,j2,i2+1));
    u4(dh4<1e-3) = nan;
    u4  = sub_chck_uv(u4,dh4,0);

% Values in the center of the grid cell	
    u34 = (u3.*dh3+u4.*dh4)./(dh3+dh4);
    u34 = u34+Vct;
    dh34= squeeze(dH(:,j2,i2));
    t34 = squeeze(T(:,j2,i2));
    s34 = squeeze(S(:,j2,i2));
    rho34= sw_dens0(s34,t34);
    dy34 = DY(j2,i2);
    hf34 = Cp*rho34.*(t34-Tref).*u34.*dh34*dy34; % W=J/s
    vf34 = u34.*dh34.*dy34;                    % m3/s - volume flux

    hflx = 0.5*(hf12+hf34);
    dZ   = 0.5*squeeze(dH(:,j1,i1)+dH(:,j2,i2));
    vflx = 0.5*(vf12+vf34);

    asgm1 = dh12*dy12;
    asgm2 = dh34*dy34;
    asgm  = 0.5*(asgm1+asgm2);
  else  % xsection
    if ny==0, error('Check norm y comp., it is 0!'); end;
    v1  = squeeze(V(:,j1,i1))*ny;
    dh1 = 0.5*(squeeze(dH(:,j1-1,i1)+dH(:,j1,i1)));
    v1(dh1<1e-3) = nan;
    v1  = sub_chck_uv(v1,dh1,0);	

    v2  = squeeze(V(:,j1+1,i1))*ny;
    dh2 = 0.5*(squeeze(dH(:,j1,i1)+dH(:,j1+1,i1)));
    v2(dh2<1e-3) = nan;
    v2  = sub_chck_uv(v2,dh2,0);	

% Values in the center of the grid cell	
    v12   = (v1.*dh1+v2.*dh2)./(dh1+dh2);
    v12   = v12+Vct;
    dh12  = squeeze(dH(:,j1,i1));	
    t12   = squeeze(T(:,j1,i1));
    s12   = squeeze(S(:,j1,i1));
    rho12 = sw_dens0(s12,t12);
    dx12  = DX(j1,i1);
    hf12  = Cp*rho12.*(t12-Tref).*v12.*dh12*dx12;
    vf12  = v12.*dh12.*dx12;                    % m3/s - volume flux

    v3  = squeeze(V(:,j2,i2))*ny;
    dh3 = 0.5*(squeeze(dH(:,j2,i2)+dH(:,j2-1,i2)));
    v3(dh3<1e-3) = nan;
    v3  = sub_chck_uv(v3,dh3,0);	

    v4  = squeeze(V(:,j2+1,i2))*ny;
    dh4 = 0.5*(squeeze(dH(:,j2,i2)+dH(:,j2+1,i2)));
    v4(dh4<1e-3) = nan;
    v4  = sub_chck_uv(v4,dh4,0);	

% Values in the center of the grid cell	
    v34 = (v3.*dh3+v4.*dh4)./(dh3+dh4);
    v34 = v34+Vct;
    dh34= squeeze(dH(:,j2,i2));	
    t34  = squeeze(T(:,j2,i2));
    s34  = squeeze(S(:,j2,i2));
    rho34= sw_dens0(s34,t34);
    dx34 = DX(j2+1,i2);
    hf34 = Cp*rho34.*(t34-Tref).*v34.*dh34*dx34;
    vf34 = v34.*dh34.*dx34;                    % m3/s - volume flux

    hflx = 0.5*(hf12+hf34);
    dZ   = 0.5*squeeze(dH(:,j1,i1)+dH(:,j2,i2));
    vflx = 0.5*(vf12+vf34);
    
    asgm1 = dh12*dx12;
    asgm2 = dh34*dx34;
    asgm  = 0.5*(asgm1+asgm2);
  end
 
%  if ig==1000,
%    fprintf('v12=%10.8d\n',v12);
%    fprintf('vf12=%10.8f, vf34=%10.8f, vflx=%10.8d, VF=%10.8d\n',...
%	  vf12(1),vf34(1),nansum(vflx),nansum(Vflx));
%  end

  Asgm(:,ig) = asgm;
  Hflx(:,ig) = hflx;
  zZ         = -cumsum(dZ);
  zZ(dZ==0)  = nan;
  ZZ(:,ig)   = zZ;
  Vflx(ig) = nansum(vflx); % depth-integrated vol flux, m3/s

end  % ig - segment

%fprintf('VFlux %4.2f Sv, HFlux %5.3d W\n',...
%	nansum(Vflx)*1e-6, nansum(nansum(Hflx)));
fprintf('VFlux=%10.8d max|Vct|=%10.8d\n',nansum(Vflx), max(max(abs(VCT))));

FSGM.Hflx_grcell_W    = Hflx;
FSGM.Vflx_Dintgr_m3_s = Vflx;
FSGM.Asgm_m2          = Asgm;
FSGM.ZZ               = ZZ;


return