% Plot surface stress emposed by sea ice on ocean
% 
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

regn = 'ARCc0.04';
expt = 22;
s_fig = 0;

%plt_fld = 'meltb'; 
%plt_fld = 'frzmlt'; 
plt_fld = 'fswdn'; % down solar flux

iyr  = 2017;
imo  = 06;
mday = 10;
btx = 'plot_seaice_fluxes.m';

pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';

ftopo = sprintf('%sdepth_ARCc0.04_17DD.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
Lmsk = HH*0;
Lmsk(HH<0)=1;


pthbin = sprintf('/nexsan/people/ddmitry/hycom/ARCc0.04_022/data/%4.4i_cice/',iyr);
fout = sprintf('%s%3.3i_cice.%i-%2.2i-%2.2i.nc',pthbin,expt,iyr,imo,mday);

fprintf('Reading %s\n',fout);


Ci = squeeze(nc_varget(fout,'aice'));
I0 = find(Ci<0.001);

% Ice heat fluxes, W/m2
% *_ai - weighted by ice area
fhocn = squeeze(nc_varget(fout,'fhocn_ai')); % heat flux to ocean, W/m2
flwup = squeeze(nc_varget(fout,'flwup_ai')); % upward longwave, 
flwdn = squeeze(nc_varget(fout,'flwdn'));    % down longwave flux
fswdn = squeeze(nc_varget(fout,'fswdn'));     % shortwave solar down
fswint= squeeze(nc_varget(fout,'fswint_ai')); % shortwave absorbed in ice interior
fswabs= squeeze(nc_varget(fout,'fswabs_ai')); % snow/ice/ocean absorbed solar flux
fsens = squeeze(nc_varget(fout,'fsens_ai')); % sensible heat flux
flat  = squeeze(nc_varget(fout,'flat_ai'));  % latent 
fswthru = squeeze(nc_varget(fout,'fswthru_ai'));
frzmlt = squeeze(nc_varget(fout,'frzmlt'));  % freeze/melt potential, >0 - new ice, <0 - ice melts

fhocn(I0)=nan;

% Ice growth: cm/day
snoice = squeeze(nc_varget(fout,'snoice')); % snow-ice formation, cm/day
congel = squeeze(nc_varget(fout,'congel')); % congelation ice growth (bottom freeze)
frazil = squeeze(nc_varget(fout,'frazil')); % frazil ice growth (1st stage of ice form)
meltt  = squeeze(nc_varget(fout,'meltt'));  % top melt
meltb  = squeeze(nc_varget(fout,'meltb'));  % basal ice melt, cm/day
meltl  = squeeze(nc_varget(fout,'meltl'));  % lateral ice melt


%[cmp,c1,c2] = sub_ice_fields(plt_fld);
switch(plt_fld),
 case('meltt');
  fld_name = 'top ice melt, cm/day';
  AA=meltt;
  c1=0;
  c2=5;

  cmp = colormap_red(400);
  cmp(1,:) = [1 1 1];

 case('meltb');
  fld_name = 'basal ice melt, cm/day';
  AA=meltb;
  c1=0;
  c2=5;

  cmp = colormap_red(400);
  cmp(1,:) = [1 1 1];

 case('frzmlt');
  fld_name = 'Freeze/melt pot. W/m^2';
  AA=frzmlt;
  c1=-400;
  c2=400;

  CMP = create_colormap3v2(400,c1,c2);
  cmp = CMP.colormap;

 case('fswdn');
  fld_name = 'ShWave Down, W/m^2';
  AA=fswdn;
  c1=0;
  c2=500;

  CMP = create_colormap6(400,c1,c2);
  cmp = CMP.colormap;

 case('flwdn');
  fld_name = 'LongWave Down, W/m^2';
  AA=flwdn;
  c1=0;
  c2=500;

  CMP = create_colormap6(400,c1,c2);
  cmp = CMP.colormap;

end;

AA(I0)=nan;

figure(1); clf;
set(gcf,'Position',[1802 637 729 694]);
contour(HH,[0 0],'k'); 
hold on
pcolor(AA); shading flat;
colormap(cmp);
caxis([c1 c2]);


axis('equal');
set(gca,'xlim',[100 3150],...
        'ylim',[500 4000],...
        'fontsize',12);

clb = colorbar;
set(clb,'Position',[0.88 0.11 0.025 0.81],...
        'Fontsize',12);

sdate = sprintf('%4.4i/%2.2i/%2.2i',iyr,imo,mday);
title(sprintf('%s, %s %s',regn,fld_name,sdate));

bottom_text(btx,'pwd',1);



