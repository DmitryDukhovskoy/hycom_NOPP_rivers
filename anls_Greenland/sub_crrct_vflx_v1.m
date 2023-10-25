function Vct = sub_crrct_vflx(GC,U,V,dH,DX,DY);
% Volume flux is not conserved
% in HYCOM mean files
% calculate net vol flux correction 
% (+Vct*dh*dz - correction for flux in 1 grid cell) 
% across the contour
% assuming it should be ~0
% disbalance distribute over
% all grid points
Hs = GC.Hbottom;
[nlr,mm,nn] = size(U);
II = GC.cntr_Iindx;
JJ = GC.cntr_Jindx;
np = length(II);
ZZ = zeros(nlr,np)*nan;
Vflx = zeros(1,np)*nan; % volume flux for checking
Asgm = zeros(1,np);
for ig=1:np-1
  i1=II(ig);
  i2=II(ig+1);
  j1=JJ(ig);
  j2=JJ(ig+1);

  if Hs(ig)>=0 & Hs(ig+1)>=0, 
    continue; 
  end;

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
  clear vf*
% collocate dH, U, V
  if xsct==0 % i2=i1
    if nx==0, error('Check norm x comp., it is 0!'); end;
    u1 = squeeze(U(:,j1,i1))*nx;
    dh1= 0.5*squeeze(dH(:,j1,i1-1)+dH(:,j1,i1));
    dy1 = DY(j1,i1);
    u1(dh1<1e-3)=nan;
    u1  = sub_chck_uv(u1);	
    vf1 = u1.*dh1.*dy1;                    % m3/s - volume flux

    u2 = squeeze(U(:,j2,i2))*nx;
    dh2= 0.5*squeeze(dH(:,j2,i2-1)+dH(:,j2,i2));
    dy2 = DY(j2,i2);
    u2(dh2<1e-3) = nan;
    u2  = sub_chck_uv(u2);	
    vf2 = u2.*dh2.*dy2;                    % m3/s - volume flux

    u3 = squeeze(U(:,j2,i2+1))*nx;
    dh3= 0.5*squeeze(dH(:,j2,i2)+dH(:,j2,i2+1));
    dy3 = DY(j2,i2+1);
    u3(dh3<1e-3) = nan;
    u3  = sub_chck_uv(u3);	
    vf3 = u3.*dh3.*dy3;                    % m3/s - volume flux

    u4 = squeeze(U(:,j1,i1+1))*nx;
    dh4= 0.5*squeeze(dH(:,j1,i1)+dH(:,j1,i1+1));
    dy4 = DY(j1,i1+1);
    u4(dh4<1e-3) = nan;
    u4  = sub_chck_uv(u4);	
    vf4 = u4.*dh4.*dy4;                    % m3/s - volume flux

    dZ   = 0.25*(dh1+dh2+dh3+dh4);
    vflx = 0.25*(vf1+vf2+vf3+vf4);
    Asgm(ig) = nansum(dZ)*0.25*(dy1+dy2+dy3+dy4);

  else  % xsection
    if ny==0, error('Check norm y comp., it is 0!'); end;
    v1  = squeeze(V(:,j1,i1))*ny;
    dh1 = 0.5*(squeeze(dH(:,j1-1,i1)+dH(:,j1,i1)));
    dx1 = DX(j1,i1);
    v1(dh1<1e-3) = nan;
    v1  = sub_chck_uv(v1);	
    vf1 = v1.*dh1.*dx1;                    % m3/s - volume flux

    v2  = squeeze(V(:,j1+1,i1))*ny;
    dh2 = 0.5*(squeeze(dH(:,j1,i1)+dH(:,j1+1,i1)));
    dx2 = DX(j1+1,i1);
    v2(dh2<1e-3) = nan;
    v2  = sub_chck_uv(v2);	
    vf2 = v2.*dh2.*dx2;                    % m3/s - volume flux

    v3  = squeeze(V(:,j2+1,i2))*ny;
    dh3 = 0.5*(squeeze(dH(:,j2,i2)+dH(:,j2+1,i2)));
    dx3 = DX(j2+1,i2);
    v3(dh3<1e-3) = nan;
    v3  = sub_chck_uv(v3);	
    vf3 = v3.*dh3.*dx3;                    % m3/s - volume flux

    v4  = squeeze(V(:,j2,  i2))*ny;
    dh4 = 0.5*(squeeze(dH(:,j2-1,i2)+dH(:,j2,i2)));
    dx4 = DX(j2,i2);
    v4(dh4<1e-3) = nan;
    v4  = sub_chck_uv(v4);	
    vf4 = v4.*dh4.*dx4;                    % m3/s - volume flux

    dZ   = 0.25*(dh1+dh2+dh3+dh4);
    vflx = 0.25*(vf1+vf2+vf3+vf4);
    Asgm(ig) = abs(nansum(dZ)*0.25*(dx1+dx2+dx3+dx4));

  end

  zZ         = -cumsum(dZ);
  zZ(dZ==0)  = nan;
  ZZ(:,ig)   = zZ;
  Vflx(ig) = nansum(vflx); % depth-integrated vol flux, m3/s

end  % ig - segment

VF = nansum(Vflx);
Hs(Hs>=0)=0;
Hz=min(ZZ);

Asct = nansum(Asgm); % lateral surface area, m2
% Correction of flux per 1 m2 of area
Vct = -VF/Asct;

%keyboard

return