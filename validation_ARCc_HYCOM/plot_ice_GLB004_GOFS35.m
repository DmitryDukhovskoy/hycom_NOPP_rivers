% Plot sea ice fields from 0.04 HYCOM-CICE
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

dnmb = datenum(2017,07,01);
DV = datevec(dnmb);



expt=216;
yr=2017;
mo=7;
dm=1;

pthmat2 = '/Net/ocean/ddmitry/HYCOM/ARCc0.08/112/data_mat/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/topo_grid/';
pthbin  = '/Net/kronos/ddmitry/hycom/GOFS3.5/';
pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc0.08/%3.3i/data_mat/',expt);
btx='plot_ice_GLB004_GOFS35.m';

% Get ARCc0.04 grid:
ftopo = sprintf('%sdepth_ARCc0.04_17DD.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
LN04 = nc_varget(ftopo,'Longitude');
LT04 = nc_varget(ftopo,'Latitude');
[mm,nn]=size(HH);




% Get ARCc grid indices
% Note that CICE grid is different from HYCOM grid
% in GOFS3.5 - over Antarctic shelf
pthycom=pthbin;
flgrd = sprintf('%sregional.grid',pthycom);
fltopo = sprintf('%sdepth_GLBc0.04_27.a',pthycom);
GRD = read_grid_bath(flgrd,fltopo);
hh=GRD.Topo;
I=find(hh>1e20);
hh=-1*hh;
hh(I)=100;


fin = sprintf('%s%3.3i_cice5_inst.%4.4i-%2.2i-%2.2i-00000.nc',...
	      pthbin,expt,yr,mo,dm);
pln = nc_varget(fin,'TLON');
I=find(pln>180);
pln(I)=pln(I)-360;
plt = nc_varget(fin,'TLAT');

% Get ARCc grid indices
ln11=LN04(1,1);
ln21=LN04(end,1);
ln12=LN04(1,end);
ln22=LN04(end,end);

lt11=LT04(1,1);
lt21=LT04(end,1);
lt12=LT04(1,end);
lt22=LT04(end,end);

dd=sqrt((pln-ln11).^2+(plt-lt11).^2);
[j11,i11]=find(dd<1e-4);
dd=sqrt((pln-ln12).^2+(plt-lt12).^2);
[j12,i12]=find(dd<1e-4);
dd=sqrt((pln-ln21).^2+(plt-lt21).^2);
[j21,i21]=find(dd<1e-4);
dd=sqrt((pln-ln22).^2+(plt-lt22).^2);
[j22,i22]=find(dd<1e-4);

IJ=[i11, j11; ...
    i12, j12; ...
    i21, j21; ...
    i22, j22];

LON = sub_Glb2Arc(pln,IJ);
LAT = sub_Glb2Arc(plt,IJ);
HH  = sub_Glb2Arc(hh,IJ);


dmm = squeeze(nc_varget(fin,'aice')); % 
hc = sub_Glb2Arc(dmm,IJ);

dmm = squeeze(nc_varget(fin,'hi'));
hi = sub_Glb2Arc(dmm,IJ);


Lmsk=HH*0;
Lmsk(HH<0)=1;


%i
stl=sprintf('GLB0.04-%3.3i, GOFS3.5, %4.4i/%2.2i/%2.2i, Ice Thickness, m',expt,yr,mo,dm);
xl1=200;
xl2=3200;
yl1=600;
yl2=4000;
c1=0;
c2=5;
fn=15;
sub_plot_hice(hi,Lmsk,xl1,xl2,yl1,yl2,c1,c2,fn);
contour(hc,[0.15 0.15],'r-');
title(stl);

bottom_text(btx,'pwd',1);

%
% Histogram of ice thicknesses
nff=2;
hbin=[0.5:10];
sub_hist_hice(hi,hbin,nff);
title(stl);

bottom_text(btx,'pwd',1,'position',[0.08 0.45 0.3 0.04]);





