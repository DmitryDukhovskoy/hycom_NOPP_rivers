% Plot statistics of the Atlantic Layer
% from GDEM climatology
% Note that this GDEM has filled in land and bottom
% Use actual topo mask 

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

%s_mat  = 0; % overridden if  s_extr=0
s_fig = 0;

%pthmat  = '/nexsan/people/takis/lydia/HYCOM/data_mat/'; 
%pthmat2 = '/Net/ocean/ddmitry/HYCOM/ARCc/data_mat/';
pthmat  = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_mat/';
%pthfig  = '/Net/ocean/ddmitry/HYCOM/ARCc/fig_GDEM4/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthin   = '/Net/data/GDEM4/';  % climatology data with no land

%fmat    = sprintf('%sADPTH_decadal.mat',pthmat); 

% Extract data for:
phi0 = 50;  % southmost latitude
zmin = -200; % depth isobath - within which perform the analysis

ftopo = sprintf('%sdepth_ARCc0.08_07.nc',pthtopo); 
HH  = nc_varget(ftopo,'Bathymetry'); 
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(HH); 

% Convert longitudes to -180 +180 rrange
% piece inserted from xsectAO.m

elon=LON; alat=LAT;

A=elon;
[my,nx]=size(elon);
for i=1:nx
for j=1:my
  long=elon(j,i);
  dmm=long;
  while (abs(dmm)>180)
    if dmm<0
      dmm=dmm+360;
    else
      dmm=dmm-360;
    end;
  end;  % while dmm
  elon(j,i)=dmm;
end;  % for j
end;  % for i
clear A;
% Grid cell spacing
[DX,DY]=sub_dx_dy(LON,LAT);


%-------------------------------
im=3;
finT=sprintf('%sptgdemv4f%2.2i.nc4',pthin,im);  % ptgdemv4f##.nc4
ZZ = -nc_varget(finT,'Depth');
TT = nc_varget(finT,'Potential_Temperature');
LNg = nc_varget(finT,'Longitude');
LT0 = nc_varget(finT,'Latitude');
iy   = max(find(LT0<=phi0))-1;
LTg  = LT0(iy+1:end);
n    = length(LNg);
m    = length(LTg);
k    = length(ZZ);
TT   = TT(:,iy+1:end,:);

% Interpolate bath into GDEM grid:
fgdm=sprintf('%sGDEM_topo.mat',pthmat);
f_intrp_Hgdem=0; % =0 - upload interpoalted, =1- interpolat, very slow

if f_intrp_Hgdem>0
  Hgdem = sub_intrp_Hgdem(LAT,LON,HH,LNg,LTg);
% Quick nearest neighbor interpolation:
  fprintf('Saving GDEM topo %s\n',fgdm);
  save(fgdm,'Hgdem');
else
  fprintf('Loading GDEM topo %s\n',fgdm);
  load(fgdm); % interpoalted for 
end

[LONg,LATg] = meshgrid(LNg,LTg);

Iarc = find(Hgdem<=-100 & LATg>70);
%Iarc = find(LATg>70); % Do over the land as well - for interpolation
nI = length(Iarc);

fprintf('Seaching for atl. water depth, GDEM, month%i ...\n',im);
Zmax = Hgdem*0;
Zt0  = Hgdem*0;
Tmax = Hgdem*0-999;

for i=1:nI
  if mod(i,1000)==0;
    fprintf('Processing, done %6.2f%%\n',i/nI*100);
  end;
  i0=Iarc(i);
  
  if (Hgdem(i0)>zmin); continue; end;
  
  tt = squeeze(TT(:,i0));
  z0 = ZZ;
  iz0 = max(find(z0>-100));
  if ~isempty(iz0)
    tt(1:iz0) = -10;
  end
  izb = max(find(z0>=Hgdem(i0)));
  tt(izb+1:end) = nan;
  
  tm = max(tt);
%    fprintf('tm=%6.4f\n',tm);
  if tm<0, continue; end;

  iz = find(tt==max(tt),1);
  Zmax(i0) = ZZ(iz);
  iz = min(find(tt>=0));

% Interpolate to find exact z
  iz1=iz-1;
  iz2=iz;
  t1 = tt(iz1);
  t2 = tt(iz2);
  z2 = ZZ(iz2);
  z1 = ZZ(iz1);
  t0 = 0;
  dtdz = (t2-t1)/(z2-z1);
  dz0 = (t0-t1)/(dtdz);
  zt0 = z1+dz0;
  Zt0(i0) = zt0;
  Tmax(i0)= tm;

end

xlim1 = 50;
xlim2 = nn-10;
ylim1 = 400;
ylim2 = mm-450;

% Interpolate GDEM into ARCc grid:
I=find(LNg>180);
LNg(I)=LNg(I)-360;
ipiv=max(find(LNg>=0));
dmm=LNg;
LNg=dmm*0;
LNg(2:n-ipiv+1,1)=dmm(ipiv+1:end);
LNg(n-ipiv+2:n+1)=dmm(1:ipiv);
dl=abs(dmm(2)-dmm(1));
LNg(1)=LNg(2)-dl;
LNg(n+2)=LNg(n+1)+dl;

% Plot 0-dgr depth
Zt1 = Zt0; 
A2=Zt0;
A2(:,2:n-ipiv+1)=Zt0(:,ipiv+1:n);
A2(:,n-ipiv+2:n+1)=Zt0(:,1:ipiv);
A2(:,1)=A2(:,n+1);
A2(:,n+2)=A2(:,2);
A2(A2==0) = nan;
A=interp2(LNg,LTg,A2,LON,LAT);
A(Hgdem>zmin)=nan;
A(A==0) = nan;
Zt0 = A;

stl = sprintf('GDEM4, Depth upper 0C m, B.Thick = 250m+/-50, Mo=%i',im);
txtb = 'plot_atlw_GDEM.m';
pfld = 'AtlZ';
nf = 1;
sub_plot_scalar(Zt0,nf,HH,xlim1,xlim2,ylim1,ylim2,LON,LAT,stl,pfld,'cmp',2);
%contour(Zt0,[-500:50:-150],'k');
%contour(Zt0,[-250 -250],'k','linewidth',1.6)
bottom_text(txtb,'pwd',1,'Fontsize',8);

if s_fig>0
  fgnm = sprintf('%sGDEM_0dgr_depth_mo%2.2i',pthfig,im);
  fprintf('Saving %s\n',fgnm);
  print('-dpng','-r300',fgnm);
end


% Plot Max T in Atl. layer:
Tmx1 = Tmax;
A2=Tmax;
A2(:,2:n-ipiv+1)=Tmax(:,ipiv+1:n);
A2(:,n-ipiv+2:n+1)=Tmax(:,1:ipiv);
A2(:,1)=A2(:,n+1);
A2(:,n+2)=A2(:,2);
A2(A2<-2) = nan;
A=interp2(LNg,LTg,A2,LON,LAT);
A(Hgdem>zmin)=nan;
%A(A==0) = nan;
Tmax = A;

stl = sprintf('GDEM4, Tmax Atl. L., B.Thick = 1C+/-0.5C, Mo=%i',im);
txtb = 'plot_atlw_GDEM.m';
pfld = 'AtlT';
nf = 2;
sub_plot_scalar(Tmax,nf,HH,xlim1,xlim2,ylim1,ylim2,LON,LAT,stl,pfld,...
                'c1',-1,'c2',4.3,'cmp',3);
%contour(Tmax,[0:0.5:2],'k');
%contour(Tmax,[1 1],'k','linewidth',1.6)
bottom_text(txtb,'pwd',1,'Fontsize',8);

if s_fig>0
  fgnm = sprintf('%sGDEM_AtlL_Tmax_mo%2.2i',pthfig,im);
  fprintf('Saving %s\n',fgnm);
  print('-dpng','-r300',fgnm);
end


% Depth of Max T below the MLD
A2=Zmax;
A2(:,2:n-ipiv+1)=Zmax(:,ipiv+1:n);
A2(:,n-ipiv+2:n+1)=Zmax(:,1:ipiv);
A2(:,1)=A2(:,n+1);
A2(:,n+2)=A2(:,2);
A2(A2==0) = nan;
A=interp2(LNg,LTg,A2,LON,LAT);
A(Hgdem>zmin)=nan;
A(A==0) = nan;
Zmax = A;

stl = sprintf('GDEM4, Depth Tmax m, Mo=%i',im);
txtb = 'plot_atlw_GDEM.m';
pfld = 'AtlTmax';
nf = 3;
sub_plot_scalar(Zmax,nf,HH,xlim1,xlim2,ylim1,ylim2,LON,LAT,stl,pfld,'cmp',2);
caxis([-600 -100]);
%contour(Zt0,[-500:50:-150],'k');
%contour(Zt0,[-250 -250],'k','linewidth',1.6)
bottom_text(txtb,'pwd',1,'Fontsize',8);






