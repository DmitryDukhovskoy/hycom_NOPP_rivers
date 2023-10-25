% Estimate dS changes for evenly distributed FW anomlay in the SPNA
% to test P. Holliday 2020 claim 6600 km3 caused 0.2 -0.3 S change
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers
startup;

close all
clear

s_fig = 0;

nbx = 5; % # of regions
regn = 'ARCc0.08';
expt = 110;  
pthfig  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/fig_NAtl/',expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
fmout   = sprintf('%sfraction_trcr_NAtl.mat',pthmat);

ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

% Grid cell spacing
[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2


xlim1 = 20;
xlim2 = nn;
ylim1 = 5;
ylim2 = 2000;

hmsk=HH;
hmsk(HH<0)=nan;

% Region of interest - double check
% with dS_FWC_timeseries_SubpolarGyre.m
%fspg='SPG_noNorth_indx.mat';  % get IGR indices of the region
%load(fspg);
%IGR(end+1,:)=IGR(1,:);

[XX,YY] = meshgrid((1:nn),(1:mm));
%INP = inpolygon(XX,YY,IGR(:,1),IGR(:,2));
%IN  = find(INP==1 & HH<0);
%INdeep = find(INP==1 & HH<-800);

% Get Subpolar gyrei
% Define Regions:
h0 = 500;
BX = sub_deep_regions(HH,LON,LAT,h0,0);
%fbx=sprintf('%sboxes_map',pthfig);
%print('-dpng','-r250',fbx);
for ib=7:7
  iBG = BX(ib).IJ;
  iBG(end+1,:) = iBG(1,:);
%  INP = inpolygon(XX,YY,iBG(:,1),iBG(:,2));
  IN = BX(ib).IN;
%  BX(ib).IN = IN;
end







% Get S fields from HYCOM
yr = 1993;
mo = 1;
dnmb = datenum(yr,mo,15);
iday = dnmb-datenum(yr,1,1)+1;
fprintf('Extracting %i/%i\n',yr,mo);
expt=110;

pthbin  = sprintf('/nexsan/archive/ARCc0.08_110/data/%4.4i/',yr);
fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);
fld='salin';
[SS,n,m,l] = read_hycom(fina,finb,fld);
%SS=squeeze(SS);
SS(SS>1e20)=nan;

[ZM,ZZ] = sub_zz_zm(fina, finb,HH);


% Need to extract monthly mean S
% for layers

  
% Volume of surplus Greenland FW in the boxes, km3
% convert to m of FW, i.e. normalize by area of region
Arg=sum(Acell(IN)); % m2
Sav=nanmean(SS(IN)); % mean S
dz0=1000;
Vanom = 6600*1e9; % FW anom, m3

% Go by layers:
clear vol 
for izz=1:41;
  zz=squeeze(ZZ(izz+1,:,:));
  zmn=nanmean(zz(IN));
  if abs(zmn)>dz0, break; end;
  dz=squeeze(abs(ZZ(izz+1,:,:)-ZZ(izz,:,:)));
  vol(izz) = nansum(Acell(IN).*dz(IN));
  DZm(izz)=nanmean(dz(IN));
end;
iz0=izz-1;
VolTot = sum(vol);  % total volume of the region 
dVdVol = Vanom./VolTot; % anom in 1 m3 of water

clear VLr Snew dS
for izz=1:iz0
  Sz = squeeze(SS(izz,:,:));
  Sav = nanmean(Sz(IN));

  VLr(izz) = nansum(dVdVol*vol(izz));
  Snew(izz) = vol(izz)*Sav/(vol(izz)+VLr(izz)); % new S
  dS (izz) = Snew(izz)-Sav;
  SAV(izz)=Sav;
end






