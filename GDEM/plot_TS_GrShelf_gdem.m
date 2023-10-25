% Plot T/S diagram for Gr Shelf
% 
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps
startup

clear
close all

disp('=======================')
disp('GDEM climatology fields');

pthin   = '/Net/data/GDEM4/';  % climatology data with no land and no bottom
pthfig  = '/Net/ocean/ddmitry/HYCOM/ARCc/fig_GDEM4/'
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = '/Net/ocean/ddmitry/HYCOM/ARCc/data_mat/';
pthmat2 = '/nexswnt/people/takis/lydia/HYCOM/data_mat/';

phi0=50;  % southmost latitude

% Get grid and bath:
%[elonA,alatA,h,my,nx,Hbth] = get_grid;   % HYCOM grid
ftopo = sprintf('%s/depth_ARCc0.08_07.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
h   = HH;
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mh,nh]=size(h); % HYCOM grid




% Winter T/S
cntr=0;
twnt=[];  % winter mean potential temperature
swnt=[];  %  winter mean S
for im=1:12
  if im>3 & im<10; continue; end;
  cntr=cntr+1;
  mnth=im;
  disp(mnth);
  iz=0;

  varnameT='Potential_Temperature'
		finT=sprintf('%sptgdemv4f%2.2i.nc4',pthin,im);  % ptgdemv4f##.nc4
		zz=nc_varget(finT,'Depth');
		elon=nc_varget(finT,'Longitude');
		alatF=nc_varget(finT,'Latitude');

		varnameS='salinity'
		finS=sprintf('%ssgdemv4f%2.2i.nc4',pthin,im);  % sgdemv4f##.nc4

%  iz=max(find(zz>=z0))-1;
  zz=-1*zz;
  iy=max(find(alatF<=phi0))-1;
  alat=alatF(iy+1:end);
  n=length(elon);
  m=length(alat);
  k=length(zz);

  disp('Reading GDEM data ...');
  PTav=nc_varget(finT,varnameT,[iz iy 0],[k m n]);
  Sav=nc_varget(finS,varnameS,[iz iy 0],[k m n]);

  if isempty(swnt)
    swnt=zeros(k,m,n);
  end
  if isempty(twnt)
    twnt=zeros(k,m,n);
  end

  twnt=twnt+PTav; % twnt is annual mean, PTav is month's contribution
  swnt=swnt+Sav;
end;   % for im
twnt=twnt./cntr;
swnt=swnt./cntr;

% GDEM bath interpolated from HYCOM ARCc topo 
fgdm=sprintf('%sGDEM_topo.mat',pthmat2);

elonM=elon;
I=elon>180;
elonM(I)=elonM(I)-360;
fprintf('Loading GDEM topo %s\n',fgdm);
load(fgdm);

% Shelf<800 m
Ish=find(Hgdem>-800 & Hgdem<0);
Inn=find(Hgdem<=-800 | Hgdem>=0);

T=twnt;
S=swnt;

T(:,Inn)=nan;
S(:,Inn)=nan;

% Specify region of interest: souther Gr shelf
lx1=1260;
lx2=1320;
ly1=35;
ly2=70;

c1=-800;
c2=0;
CMP=create_colormap6(400,c1,c2);
cmp=CMP.colormap;
cnt=CMP.intervals;

T=T(:,ly1:ly2,lx1:lx2);
S=S(:,ly1:ly2,lx1:lx2);
nz=length(zz);

%
% Create dens
tt=[-1:0.1:8];
ss=[32:0.1:37];
[Sb,Tb]=meshgrid(ss,tt);
rhoW=sw_dens0(Sb,Tb)-1000;


figure(1); clf;
hold on;

for iz=1:nz
  z0=zz(iz);
  iclr=min(find(cnt>=z0));
  clr=cmp(iclr,:);

  tmm=squeeze(T(iz,:,:));
  smm=squeeze(S(iz,:,:));
  I0=find(~isnan(tmm));
  tmm=tmm(I0);
  smm=smm(I0);
 
  plot(tmm,smm,'.','Color',clr,'Markersize',10);
end

contour(Tb,Sb,rhoW,[25.6:0.2:30],'Color',[0.6 0.6 0.6]);
contour(Tb,Sb,rhoW,[28 28],'Color',[0.3 0.3 0.3]);
contour(Tb,Sb,rhoW,[26 26],'Color',[0. 0.4 0.6]);
set(gca,'xlim',[0 6.2],...
        'ylim',[33 36],...
        'tickdir','out',...
        'xtick',[0:0.5:8],...
        'ytick',[32:0.5:36],...
        'Fontsize',14);

clb=colorbar;
colormap(cmp);
caxis([c1 c2]);
set(clb,'position', [0.91 0.26 0.015 0.5],...
        'Fontsize',14);

title('T-S SW Gr Shelf, Winter, GDEMv4');








