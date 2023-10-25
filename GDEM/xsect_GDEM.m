% Calculate mixing layer depth in the Arctic ocean 
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

pthin   = '/Net/data/GDEM4/';  % climatology data with no land and no bottom
pthfig  = '/Net/ocean/ddmitry/HYCOM/ARCc/fig_GDEM4/'
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = '/Net/ocean/ddmitry/HYCOM/ARCc/data_mat/';
pthmat2 = '/nexsan/people/takis/lydia/HYCOM/data_mat/';

% comment out the unnecessary stuff
%var='T';
%target=0; cmin=-2;cmax=8;
%target=-250; cmin=-2; cmax=8;
%target=-1000; cmin=-1; cmax=0;
%units='C';
%
%var='S';
%target=0; cmin=28;cmax=36;
%target=-250; cmin=33; cmax=35;
%target=-1000; cmin=34.85; cmax=34.95;
%units='ppt';

   

% Flags
sfig  = 1;  % 1 - save figures

% Extract data for:
phi0=50;  % southmost latitude
%year=2006;

% Get grid and bath:
%[elonA,alatA,h,my,nx,Hbth] = get_grid;   % HYCOM grid
ftopo = sprintf('%s/depth_ARCc0.08_07.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
h   = HH;
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mh,nh]=size(h); % HYCOM grid



cntr=0;
ptan=[];  % annual mean potential temperature
san=[];  % annual mean salinity
for im=1:1
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

  if isempty(san)
    san=zeros(k,m,n);
  end
  if isempty(ptan)
    ptan=zeros(k,m,n);
  end

  ptan=ptan+PTav; % ptan is annual mean, PTav is month's contribution
  san=san+Sav; 
end;   % for im
ptan=ptan./cntr;
san=san./cntr;

% Interpolate bath into GDEM grid:
fgdm=sprintf('%sGDEM_topo.mat',pthmat2);

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

% set layers deeper than bottom to nan
%    for kk=1:k
%    % set layers over land to nan;
%        ptan(kk,Hgdem>=0)=nan;
%        san(kk,Hgdem>=0)=nan;
%    % if Hgdem(i,j)>zz(kk) set ptan(kk,i,j)=nan;
%        ptan(kk,Hgdem>zz(kk))=nan;
%  	san(kk,Hgdem>zz(kk))=nan;
%    end
%

%keyboard

 x1lon=-150; x2lon=70;
 x3lon=-100; x4lon=80;

 x1lat=alat;x2lat=fliplr(alat);
 x3lat=alat;x4lat=fliplr(alat);

 x12lat=[x1lat;x2lat];
 x34lat=[x2lat;x3lat];

 x1san=squeeze(san(:,:,find(elonM==x1lon))); x2san=fliplr(squeeze(san(:,:,find(elonM==x2lon))));
 x3san=squeeze(san(:,:,find(elonM==x3lon))); x4san=fliplr(squeeze(san(:,:,find(elonM==x4lon))));

 x12san=[x1san,x2san];x34san=[x3san,x4san];

 x1ptan=squeeze(ptan(:,:,find(elonM==x1lon))); x2ptan=fliplr(squeeze(ptan(:,:,find(elonM==x2lon))));
 x3ptan=squeeze(ptan(:,:,find(elonM==x3lon))); x4ptan=fliplr(squeeze(ptan(:,:,find(elonM==x4lon))));

 x12ptan=[x1ptan,x2ptan];x34ptan=[x3ptan,x4ptan];

 close all;
 figure(1)
 pcolor(1:size(x12lat),zz,x12ptan); shading flat; colorbar
 figure(2)
 pcolor(1:size(x12lat),zz,x12san); shading flat; colorbar
 figure(3)
 pcolor(1:size(x34lat),zz,x34ptan); shading flat; colorbar
 figure(4)
 pcolor(1:size(x34lat),zz,x34san); shading flat; colorbar


 













