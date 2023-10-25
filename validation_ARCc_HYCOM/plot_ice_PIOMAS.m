addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close

s_mat_indx = 0; % =1 - derive indices/weights for interpolation, save mat
              % =0 - load saved mat with weights
s_mat   = 0;  % =1 - save PIOMAS fields interpolated into HYCOM 

PTH.topo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
PTH.data = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/PIOMAS_ice_data/';
PTH.rest = '/nexsan/people/ddmitry/hycom/ARCc0.08/103/ice_restart/';
PTH.mat  = PTH.rest;

btx = 'plot_ice_PIOMAS.m'; 


fgrds = sprintf('%sgrid.dat',PTH.data); % grid for scalar fields
nxl=360;
nyl=120;
nxy=nxl*nyl;

% read lon/lat scalar fields
dmm = load(fgrds);
[a1,a2]=size(dmm);
nrw=nxy/a2;

LONp = dmm(1:nrw,:);
LONp = reshape(LONp',[nxy,1]);
LONp = reshape(LONp,[nxl,nyl])';
[mp,np]=size(LONp);

LATp = dmm(nrw+1:end,:);
LATp = reshape(LATp',[nxy,1]);
LATp = reshape(LATp,[nxl,nyl])';

% HYCOM:
fltopo=sprintf('%sdepth_ARCc0.08_11.nc',PTH.topo);
HH   = nc_varget(fltopo,'Bathymetry');
alat = nc_varget(fltopo,'Latitude');
elon = nc_varget(fltopo,'Longitude');
LATh  = alat;
LONh  = elon;
[m,n]= size(HH);
[mm,nn]= size(HH);

Lmsk=HH*0;
Lmsk(HH<0)=1;

% Grid cell spacing
[DX,DY]=sub_dx_dy(elon,alat);
Acell=DX.*DY; % Grid cell area, m2


H2P = sub_intrp_indx_piomas2hycom(s_mat_indx);

% Read PIOMAS:
% For scalar data use LONp/LATp


yr=2017;
imo=7;
iday=1;

%fld='area';
fld='aiday';
imoday=imo*100+iday;
APIOMAS = sub_read_piomas(fld,yr,imoday);
FP = APIOMAS.Field;
% Interpolate PIOMAS to HYCOM
FF = HH*0;
WW = H2P.Weights;
IP = H2P.PIOMAS_LinIndx;
IH = H2P.HYCOM_IceLinIndx;
dmm = FP(IP(:,1)).*WW(:,1)+...
      FP(IP(:,2)).*WW(:,2)+...
      FP(IP(:,3)).*WW(:,3)+...
      FP(IP(:,4)).*WW(:,4)+...
      FP(IP(:,5)).*WW(:,5);
FF(IH)=dmm;
FF(FF==0)=nan;
hcnc=FF;




%fld='ithkn'; % monthly mean
%APIOMAS = sub_read_piomas(fld,yr,imo);
%FP = APIOMAS.Field;
fld='hiday';  % daily mean
imoday=imo*100+iday;
APIOMAS = sub_read_piomas(fld,yr,imoday);
FP = APIOMAS.Field;

% Interpolate PIOMAS to HYCOM
FF = HH*0;
WW = H2P.Weights;
IP = H2P.PIOMAS_LinIndx;
IH = H2P.HYCOM_IceLinIndx;
dmm = FP(IP(:,1)).*WW(:,1)+...
      FP(IP(:,2)).*WW(:,2)+...
      FP(IP(:,3)).*WW(:,3)+...
      FP(IP(:,4)).*WW(:,4)+...
      FP(IP(:,5)).*WW(:,5);
FF(IH)=dmm;
%FF(FF==0)=nan;
hice=FF;

hice(Lmsk==0)=nan;

stl=sprintf('PIOMAS Reanls, %4.4i/%2.2i/%2.2i, Ice Thickness, m',yr,imo,iday);
xl1=100;
xl2=1600;
yl1=300;
yl2=2000;
c1=0;
c2=5;
fn=15;
sub_plot_hice(hice,Lmsk,xl1,xl2,yl1,yl2,c1,c2,fn);
contour(hcnc,[0.15 0.15],'r-');
title(stl);

bottom_text(btx,'pwd',1);

%
% Histogram of ice thicknesses
sttl=sprintf('PIOMAS, Hice, %4.4i/%2.2i/%2.2i, Ice Thickness, m',yr,imo,iday);
nff=2;
hbin=[0.5:10];
sub_hist_hice(hice,hbin,nff);
title(sttl);

bottom_text(btx,'pwd',1,'position',[0.08 0.45 0.3 0.04]);





