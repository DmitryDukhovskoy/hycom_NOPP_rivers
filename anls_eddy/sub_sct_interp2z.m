function UTSZ = sub_sct_interp2z(SCT,U,V,T,S,ZZf,ZZh,HH,dH,DX,DY);
%  NOTE : bin-averaging for z interpolation needs further checking
%  not recommended to use this code
 
%
% Interpolate normal to a vertical section 
% (2D array)
% SCT - structured aray with info about the section
%
% onto z-fixed depths (ZZf)
% ZZh - HYCOM 3D array of layer interface depths
% ZZf - 1D array of depths fixed depths
%
% Approach:
% Averaging:
%            v1(i,j+1)        v3(i+1,j+1)
%      |--------|-------|-------|---------|
%      |                |                 |
%      |                |                 |
%      |     V12        |       V34       | Flx = Cp*rho*(T-Tref)V1*dH1*dX1/2+
%      -        *=======X========*        |       Cp*rho*(T-Tref)V2*dH2*dX2/2+
%u(i,j)|     T,S,dH     |                 | where V1 and V2 are flux-averaged
%      |                |                 | v(i,j+1)&v(i,j) 
%      |                |                 |  and v(i+1,j+1)&v(i+1,j)
%      |                |                 |
%      |--------|-------|-------|---------|
%      |      v2(i,j)         v4(i+1,j)
%      |   
%
%  Average and collocate all variables into "X" point
%  Then interpolate onto regular grid
fprintf('Interpolating U,T,S Greenl. contour -> z\n');

hg = 1e20;
[nlr,mm,nn] = size(U);
nZ = length(ZZf);
mZ = length(ZZf)-1;
%keyboard

nsct=length(SCT);
for js=1:nsct
%  fprintf('Segment %i\n',js);
  II = SCT(js).I;
  JJ = SCT(js).J;
  Hs = SCT(js).Hbottom;

  np = length(II);
  Uh = zeros(nlr,np)+hg;
  Th = zeros(nlr,np)+hg;
  Sh = zeros(nlr,np)+hg;
  UUh= zeros(nlr,np)+hg;
  VVh= zeros(nlr,np)+hg;
  Zh = zeros(nlr+1,np)+hg;
  Uz = zeros(mZ,np)+hg;
  Tz = zeros(mZ,np)+hg;
  Sz = zeros(mZ,np)+hg;
  UUz= zeros(mZ,np)+hg;
  VVz= zeros(mZ,np)+hg;
  

  for ig=1:np-1  % small segments
    i1=II(ig);
    i2=II(ig+1);
    j1=JJ(ig);
    j2=JJ(ig+1);

    clear u12 u34 v12 v34 

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

    nrm = SCT(js).Nrm(ig,:);
    nx  = nrm(1);
    ny  = nrm(2);
    dZ = 0;
    clear hf* vf*

    if xsct==0 % i2=i1 Y-section, use u components
      if nx==0, error('Check norm x comp., it is 0!'); end;
  % U Fluxes in grid cell (j1,i1):	
      u1 = squeeze(U(:,j1,i1))*nx;
      dh1= 0.5*squeeze(dH(:,j1,i1-1)+dH(:,j1,i1));
      u1(dh1<1e-3)=nan;
      u1  = sub_chck_uv(u1,dh1,0);
      u1(isnan(u1))=0;

      u2 = squeeze(U(:,j1,i1+1))*nx;
      dh2= 0.5*squeeze(dH(:,j1,i1)+dH(:,j1,i1+1));
      u2(dh2<1e-3) = nan;
      u2  = sub_chck_uv(u2,dh2,0);
      u2(isnan(u2))=0;

      I=find(dh1==0 & dh2==0);
      u1(I)=nan;
      u2(I)=nan;

  % U Fluxes in grid cell (j2,i2):	
      u3  = squeeze(U(:,j2,i2))*nx;
      dh3 = 0.5*squeeze(dH(:,j2,i2-1)+dH(:,j2,i2));
      u3(dh3<1e-3) = nan;
      u3  = sub_chck_uv(u3,dh3,0);
      u3(isnan(u3))=0;

      u4  = squeeze(U(:,j2,i2+1))*nx;
      dh4 = 0.5*squeeze(dH(:,j2,i2)+dH(:,j2,i2+1));
      u4(dh4<1e-3) = nan;
      u4  = sub_chck_uv(u4,dh4,0);
      u4(isnan(u4))=0;

      I=find(dh3==0 & dh4==0);
      u3(I)=nan;
      u4(I)=nan;

    else % X-section, V components
      u1 = squeeze(V(:,j1,i1))*ny;
      dh1= 0.5*squeeze(dH(:,j1-1,i1)+dH(:,j1,i1));
      u1(dh1<1e-3)=nan;
      u1  = sub_chck_uv(u1,dh1,0);
      u1(isnan(u1))=0;

      u2  = squeeze(V(:,j1+1,i1))*ny;
      dh2 = 0.5*(squeeze(dH(:,j1,i1)+dH(:,j1+1,i1)));

      u2(dh2<1e-3) = nan;
      u2  = sub_chck_uv(u2,dh2,0);
      u2(isnan(u2))=0;

      I=find(dh1==0 & dh2==0);
      u1(I)=nan;
      u2(I)=nan;

      u3  = squeeze(V(:,j2,i2))*ny;
      dh3 = 0.5*(squeeze(dH(:,j2,i2)+dH(:,j2-1,i2)));
      u3(dh3<1e-3) = nan;
      u3  = sub_chck_uv(u3,dh3,0);
      u3(isnan(u3))=0;

      u4  = squeeze(V(:,j2+1,i2))*ny;
      dh4 = 0.5*(squeeze(dH(:,j2,i2)+dH(:,j2+1,i2)));
      u4(dh4<1e-3) = nan;
      u4  = sub_chck_uv(u4,dh4,0);
      u4(isnan(u4))=0;

      I=find(dh3==0 & dh4==0);
      u3(I)=nan;
      u4(I)=nan;

    end

  % Values in the center of the grid cell	
    u12 = (u1.*dh1+u2.*dh2)./(dh1+dh2);
    dh12= squeeze(dH(:,j1,i1));
    t12 = squeeze(T(:,j1,i1));
    s12 = squeeze(S(:,j1,i1));
    dy12 = DY(j1,i1);

  % Values in the center of the grid cell	
    u34 = (u3.*dh3+u4.*dh4)./(dh3+dh4);
    dh34= squeeze(dH(:,j2,i2));
    t34 = squeeze(T(:,j2,i2));
    s34 = squeeze(S(:,j2,i2));
    dy34 = DY(j2,i2);

  % Estimate all variables in the middle of the segment
    u12(isnan(u12))=0;
    u34(isnan(u34))=0;
    I=find(dh12==0 & dh34==0);
  %  u12(I)=nan;
  %  u34(I)=nan;
  %  t12(I)=nan;
  %  t34(I)=nan;
  %  s12(I)=nan;
  %  s34(I)=nan;
    dh0 = 0.5*(dh12+dh34);
    u0 = (u12.*dh12+u34.*dh34)./(dh12+dh34);
    t0 = (t12.*dh12+t34.*dh34)./(dh12+dh34);
    s0 = (s12.*dh12+s34.*dh34)./(dh12+dh34);

    hbtm = 0.5*(Hs(ig)+Hs(ig+1));
    zh12 = squeeze(ZZh(:,j1,i1));
    zh34 = squeeze(ZZh(:,j2,i2));
    zh00  = 0.5*(zh12+zh34);
  % correct L. thickness to match
  % total bottom depth
    zz1=0;
    zz2=0;
    for kll=1:nlr-1
      dz1=abs(dh0(kll));
      dz2=abs(dh0(kll+1));
      zz1=zz1-dz1;
      zz2=zz1-dz2;
      if abs(zz2)>abs(hbtm),
	dz2=hbtm-zz1;
	dh0(kll+1:end)=0;
	dh0(kll+1)=abs(dz2);
	break;
      end
    end

    zh0=zh00*0;
    for kll=1:nlr
      zh0(kll+1,1)=zh0(kll)-dh0(kll);
    end

  % Check:
    errZ = max(abs(zh00-zh0));
    if errZ>0.1
      fprintf('sub_interp2z: Error in mid-point depth array errZ=%8.6f\n',errZ);
      fprintf('Check local depth, layer thicknesses: dh0,zh0\n');
      fprintf('=======   STOPPING     ======= \n');
      keyboard
    end

% Extract U,V components at the segment vertices
% for plotting vectors along the sections
    uu=squeeze(U(:,j1,i1));
    vv=squeeze(V(:,j1,i1));
    
    fchck=0;
    u0i = sub_binav2z_1D(hbtm,zh0,ZZf,u0,fchck);
    t0i = sub_binav2z_1D(hbtm,zh0,ZZf,t0,fchck);
    s0i = sub_binav2z_1D(hbtm,zh0,ZZf,s0,fchck);
    uui = sub_binav2z_1D(hbtm,zh0,ZZf,uu,fchck);
    vvi = sub_binav2z_1D(hbtm,zh0,ZZf,vv,fchck);

  %  if Hs(ig)<=-1000,
  %    keyboard
  %  end

    Uz(:,ig)=u0i;
    Tz(:,ig)=t0i;
    Sz(:,ig)=s0i;
    UUz(:,ig)=uui;
    VVz(:,ig)=vvi;

    Uh(:,ig)=u0;
    Th(:,ig)=t0;
    Sh(:,ig)=s0;
    Zh(:,ig)=zh0;
    UUh(:,ig)=uu;
    VVh(:,ig)=vv;

  end  % segments

  Uz(Uz>0.1*hg)=nan;
  Tz(Tz>0.1*hg)=nan;
  Sz(Sz>0.1*hg)=nan;
  UUz(UUz>0.1*hg)=nan;
  VVz(VVz>0.1*hg)=nan;
  
  Uh(Uh>0.1*hg)=nan;
  Th(Th>0.1*hg)=nan;
  Sh(Sh>0.1*hg)=nan;
  Zh(Zh>0.1*hg)=nan;
  UUh(UUh>0.1*hg)=nan;
  VVh(VVh>0.1*hg)=nan;

  UTSZ(js).NormalU=Uz;   % normal component in the middle of gr cells
  UTSZ(js).Temp=Tz;      % T in the middle
  UTSZ(js).Salin=Sz;     % S in the middle
  UTSZ(js).Ucomp=UUz;    % U comp at gr cell (i,j) not collocated
  UTSZ(js).Vcomp=VVz;    % V comp not collocated
  
  f_chck=0;
  if f_chck==1 & js==3
    X = [1:np];
    [XX,dmm]=meshgrid(X,[1:nlr]);
    zh = Zh(1:end-1,:);
    
%    ah = Uh;
%    az = Uz;
%  c1=-0.5;
%  c2=0.5;

  %  ah = Th;
  %  az = Tz;
  %  c1=-1;
  %  c2=4;
  
    ah = Sh;
    az = Sz;
    c1=33;
    c2=35;

    figure(11); clf;
    pcolor(XX,zh,ah); shading flat;
    colorbar;
    caxis([c1 c2]);
    set(gca,'ylim',[-1200 0]);
    title('Original Field, HYCOM grid');


    figure(12); clf
    pcolor(X,ZZf(1:mZ),az); shading flat;
    colorbar;
    caxis([c1 c2]);
    set(gca,'ylim',[-1200 0]);
    title('Interpolated Z');

     keyboard

  end

end; % segments





return
