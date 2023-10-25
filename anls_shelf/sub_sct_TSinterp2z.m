function TSZ = sub_sct_TSinterp2z(SCT,T,S,ZZf,ZZh,HH,dH,DX,DY);

% Interpolate T&S along the sections
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
fprintf('Interpolating T,S vertical section  -> z levels\n');

hg = 1e20;
[nlr,mm,nn] = size(T);
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
%  Uh = zeros(nlr,np)+hg;
  Th = zeros(nlr,np)+hg;
  Sh = zeros(nlr,np)+hg;
%  UUh= zeros(nlr,np)+hg;
%  VVh= zeros(nlr,np)+hg;
  Zh = zeros(nlr+1,np)+hg;
%  Uz = zeros(mZ,np)+hg;
  Tz = zeros(mZ,np)+hg;
  Sz = zeros(mZ,np)+hg;
%  UUz= zeros(mZ,np)+hg;
%  VVz= zeros(mZ,np)+hg;
  

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


  % Values in the center of the grid cell	
    dh12= squeeze(dH(:,j1,i1));
    t12 = squeeze(T(:,j1,i1));
    s12 = squeeze(S(:,j1,i1));
    dy12 = DY(j1,i1);

  % Values in the center of the grid cell	
    dh34= squeeze(dH(:,j2,i2));
    t34 = squeeze(T(:,j2,i2));
    s34 = squeeze(S(:,j2,i2));
    dy34 = DY(j2,i2);

  % Estimate all variables in the middle of the segment
    I=find(dh12==0 & dh34==0);
    dh0 = 0.5*(dh12+dh34);
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

    
    fchck=0;
    t0i = sub_binav2z_1D(hbtm,zh0,ZZf,t0,fchck);
    s0i = sub_binav2z_1D(hbtm,zh0,ZZf,s0,fchck);

  %  if Hs(ig)<=-1000,
  %    keyboard
  %  end

    Tz(:,ig)=t0i;
    Sz(:,ig)=s0i;

    Th(:,ig)=t0;
    Sh(:,ig)=s0;
    Zh(:,ig)=zh0;

  end  % segments

  Tz(Tz>0.1*hg)=nan;
  Sz(Sz>0.1*hg)=nan;
  
  Th(Th>0.1*hg)=nan;
  Sh(Sh>0.1*hg)=nan;
  Zh(Zh>0.1*hg)=nan;

  TSZ(js).Temp=Tz;
  TSZ(js).Salin=Sz;
  
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