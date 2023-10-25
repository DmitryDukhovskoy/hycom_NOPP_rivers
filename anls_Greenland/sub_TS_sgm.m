function FSGM = sub_TS_sgm(GC,S,T,dH,DX,DY,HH);
% Calculate vol and heat fluxes across segments
% of the contour

[nlr,mm,nn] = size(S);

II = GC.cntr_Iindx;
JJ = GC.cntr_Jindx;
Hs = GC.Hbottom;

np   = length(II);
Tsgm = zeros(nlr,np)*nan;
Ssgm = zeros(nlr,np)*nan;
%Rsgm = zeros(nlr,np)*nan;
ZZ   = zeros(nlr,np)*nan;
DZ   = zeros(nlr,np)*nan;
for ig=1:np
  i1=II(ig);
  j1=JJ(ig);

  if HH(j1,i1)>=0, 
    continue; 
  end;

  if Hs(ig)>0, continue; end;

  dZ  = 0;
  t1  = squeeze(T(:,j1,i1));
  s1  = squeeze(S(:,j1,i1));
%  rho1= sw_dens0(s1,t1);
  dZ  = squeeze(dH(:,j1,i1));

  Tsgm(:,ig) = t1;
  Ssgm(:,ig) = s1;
%  Rsgm(:,ig) = rho1;
  zZ         = -cumsum(dZ);
  zZ(dZ==0)  = nan;
  ZZ(:,ig)   = zZ;
  DZ(:,ig) = dZ;

end  % ig - segment

%fprintf('VFlux %4.2f Sv, HFlux %5.3d W\n',...
%	nansum(Vflx)*1e-6, nansum(nansum(Hflx)));
%fprintf(':::  Min T= %6.4f, Max T=%6.4f\n',min(min(Tsgm)),max(max(Tsgm)));
%fprintf(':::  Min S= %6.4f, Max S=%6.4f\n',min(min(Ssgm)),max(max(Ssgm)));

FSGM.Tsgm    = Tsgm;
FSGM.Ssgm    = Ssgm;
FSGM.ZZ      = ZZ;
FSGM.DZ      = DZ;


return