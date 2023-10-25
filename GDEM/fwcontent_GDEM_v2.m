% Calculate freshwater content in the Arctic Ocean
%  GDEM climatology - annual
% Note that FWC is integrated from the isohaline = Sref up to the surface
% See Haine et al. 2014 "Arctic FW export ...", Global and planetary change

addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps
startup

clear 
close all

disp('=======================')
disp('GDEM climatology fields');

%pthin='/gfs3/hycom_dataserver/ftp/datasets/ieee/clim/GDEM3/';  % climatology data with no land
pthin   = '/Net/data/GDEM3/';  % climatology data with no land
pthfig  = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.72/GDEM/fig/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = '/Net/tholia/ddmitry/hycom/ARCc0.08/data_mat/';

% Flags
sfig  = 0;  % 1 - save figures

% Extract data for:
phi0=50;  % southmost latitude
%year=2006;

% Get grid and bath:
%[elonA,alatA,h,my,nx,Hbth] = get_grid;   % HYCOM grid
ftopo = sprintf('%s/depth_ARCc0.08_09.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
h   = HH;
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mh,nh]=size(h); % HYCOM grid



cntr=0;
san=[];
for im=1:12
  cntr=cntr+1;
  mnth=im;
  disp(mnth);
  fin=sprintf('%ss_m%2.2i.nc4',pthin,im);
  zz=nc_varget(fin,'depth');
  zz=-1*zz;
%  iz=max(find(zz>=z0))-1;
  elon=nc_varget(fin,'longitude');
  alatF=nc_varget(fin,'latitude');
  iy=max(find(alatF<=phi0))-1;
  alat=alatF(iy+1:end);
  n=length(elon);
  m=length(alat);
  k=length(zz);
  
  disp('Reading GDEM data ...');
  iz=0;
  Sav=nc_varget(fin,'Salinity',[iz iy 0],[k m n]);
  if isempty(san)
    san=zeros(k,m,n);
  end

  san=san+Sav;
end;   % for im
san=san./cntr;

% Interpolate bath into GDEM grid:
fgdm=sprintf('%sGDEM_topo.mat',pthmat);

elonM=elon;
I=elon>180;
elonM(I)=elonM(I)-360;

f_intrp_Hgdem=0; % =0 - upload interpoalted, =1- interpolat, very slow

if f_intrp_Hgdem>0
  Hgdem = sub_intrp_Hgdem(LAT,LON,HH,elon,alat);
% Quick nearest neighbor interpolation:
  fprintf('Saving GDEM topo %s\n',fgdm);
  save(fgdm,'Hgdem');
else
  fprintf('Loading GDEM topo %s\n',fgdm);
  load(fgdm);
end
  




Sref=34.8;
zw(1)=0;
a1=length(zz);
for kk=1:a1-1
  dmm=0.5*abs(zz(kk+1)-zz(kk));
  zw(kk+1)=zz(kk)-dmm;
end;
zw(a1+1)=zz(a1)-dmm;

% Integrate over depth
disp('Calculating FWC ...');
zmin=-500;
%zint=-500;  % depth should be - depth of Sref. isohaline.
%dmm=abs(zz-zint);
%kz=find(dmm==min(dmm));
kz=[];
Fwc=zeros(m,n)*nan;
for i=1:n
  for j=1:m
    if Hgdem(j,i)>zmin, 
      continue; 
    end; % 
    dmm=san(1,j,i);
    
% Integrate from the depth of Sref
% or from the bottom, if shallower
    Sz=squeeze(san(:,j,i));
    kz=max(find(Sz<=Sref));
    if isempty(kz) % no freshwater
      Fwc(j,i)=0;
      continue;
    end
    
%    keyboard
  
    if kz>0,
      asum=0;
      for k=1:kz
	dz=abs(zw(k+1)-zw(k));
  %      if (san(k,j,i)<=Sref)
	asum=asum+dz*((Sref-san(k,j,i))/Sref);
  %      end;
      end;
      dmm=asum;
      Fwc(j,i)=dmm;
    else
      Fwc(j,i)=0;
    end;  % if
  %  Fwc(j,i)=dmm;
  
  end;
end;
%I=find(Fwc<1);
%Fwc(I)=nan;
% 
% Interpolate GDEM into ARCc grid:
    I=find(elon>180);
    elon(I)=elon(I)-360;
    ipiv=max(find(elon>=0));
    dmm=elon;
    elon=dmm*0;
    elon(2:n-ipiv+1,1)=dmm(ipiv+1:end);
    elon(n-ipiv+2:n+1)=dmm(1:ipiv);
    dl=abs(dmm(2)-dmm(1));
    elon(1)=elon(2)-dl;
    elon(n+2)=elon(n+1)+dl;

  A2=Fwc;
  A2(:,2:n-ipiv+1)=Fwc(:,ipiv+1:n);
  A2(:,n-ipiv+2:n+1)=Fwc(:,1:ipiv);
  A2(:,1)=A2(:,n+1);
  A2(:,n+2)=A2(:,2);
  A=interp2(elon,alat,A2,LON,LAT);
  A(HH>zmin)=nan;
  A(:,1:300)=nan;
  A(:,1300:end)=nan;
  A(1:350,:)=nan;
  A(2000:end,:)=nan;
%  MLD=A;

%keyboard

% Colormap:
c1=-10;
c2=20;
cl1=colormap_cold(20);
cl2=flipud(colormap_red(10));
cmp=[cl2;cl1];
nint=length(cmp);
cnt=[c1:(c2-c1)/nint:c2];


hmsk=HH;
hmsk(HH<0)=1;
hmsk(HH>0)=0;


cmp0=[0,0,0;
      0.7, 0.7, 0.7];
figure(1); clf;
axes('Position',[0.08 0.2 0.87 0.7]);
pcolor(hmsk); shading flat;
colormap(cmp0);
axis('equal');
set(gca,'xlim',[50 1550],...
	'ylim',[400 2000],...
	'xtick',[],...
	'ytick',[]);
freezeColors;

stt=sprintf('FW Content (m) from GDEM, Sref=%5.2f, z0=%5.1f',Sref,zmin);
title(stt,'Fontsize',14);



ps=get(gca,'Position');

axes('Position',ps);
pcolor(A); shading flat;
caxis([c1 c2]);
colormap(cmp);
hold on;

nf=1;
dlmb=45;
dphi=10;
clr=[0.8 0.8 0.8];
plot_gridlines(dlmb,dphi,nf,clr,LON,LAT);

axis('equal');
set(gca,'xlim',[50 1550],...
	'ylim',[400 2000],...
	'xtick',[],...
	'ytick',[]);
set(gca,'visible','off');


%axes('Position',[0.08 0.04 0.87 0.1]);
hght=[];
lngth=[];
mint=2;
mbx=mint;
fsz=14;
bxc='w';
%    posc=[0.23 0.04 0.6 0.05];
posc=[0.23 0.14 0.6 0.05];
nsint2=mint;
aend=1;
[az,axc]  = colorbar_horiz (cmp,cnt,hght,lngth,mint,fsz,bxc,posc,mbx,aend);

%set(gca,'visible','off');
txtb='/hycom_arc72/GDEM_analysis/fwcontent_GDEM.m';
bottom_text(txtb);

if sfig==1
  fgnm=sprintf('%sGDEM_FWcont_annual_z0-%im',pthfig,abs(zint));
  fprintf('Saving %s\n',fgnm);
  print('-dpng',fgnm);
end














