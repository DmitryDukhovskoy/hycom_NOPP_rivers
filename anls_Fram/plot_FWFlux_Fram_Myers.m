addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

YR1=2002;
YR2=2016;

pthdat='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/data_Myers/FRAM/';

btx = 'calc_FWFlux_Fram_Myres.m';

iSp=199; % index of Spitsbergen coast
for iyr=YR1:YR2
  dJ1=datenum(iyr,1,1);
  for iday=5:5:365
    dnmb=dJ1+iday-1;
    dv=datevec(dnmb);
    mo=dv(2);
    mday=dv(3);
    fnm=sprintf('%sANHA12-EXH006_y%4.4im%2.2id%2.2i_FRAMclip_gridV.nc',...
	      pthdat,iyr,mo,mday);
    
    if ~exist('ZZ','var');
      ZZ=-nc_varget(fnm,'depthv');
      xlon=nc_varget(fnm,'nav_lon');
      xlat=nc_varget(fnm,'nav_lat');
      nz=length(ZZ);
      nx=length(xlon);
    end
    
    V=squeeze(nc_varget(fnm,'vomecrty'));
    
    pcolor(xlon,ZZ,V); shading flat;
    
    
    