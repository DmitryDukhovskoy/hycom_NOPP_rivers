% Calculate volume of CICE output - FWC in solid phase
% for the same regions as oceanic FWC 
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

expt=110;
s_mat = 1;

yr=1993;
mo=1;
dm=2;

%pthbin  = '/Net/kronos/ddmitry/hycom/ARCc0.08/tmp_files/';
%pthbin  = '/nexsan/people/ddmitry/hycom/ARCc0.08/110/cice/';
pthbin = sprintf('/nexsan/archive/ARCc0.08_%s/data/%i_cice/',expt,yr);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/data_mat/';

fmat    = sprintf('%scice_vol_arc08_%3.3i.mat',pthmat,expt);

ftopo = sprintf('%sdepth_ARCc0.08_11.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(HH);

fin = sprintf('%s%3.3i_cice_inst.%4.4i-%2.2i-%2.2i-00000.nc',...
	      pthbin,expt,yr,mo,dm);

% Grid cell spacing
[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2

% Define indices of the select regions:
% See anls_Greenland sub_define_boxes
FIBX = sub_define_boxesAO_Natl(HH,LON,LAT,0);
nbx = length(FIBX);

% Monthly values
YRPLT=[];
cc=0;
for yr=1993:2015
  for mo=1:12
    for mday=5:10:30
      cc=cc+1;
      YRPLT(cc,1)=yr;
      YRPLT(cc,2)=mo;
      YRPLT(cc,3)=mday;
    end
  end
end
nrec=cc;

%Start time loop:
cc  = 0;
cyr = 0;
ndm = 3; % #days per mont
cday= 0;
Fwc1Yr = HH*0;
Fwc2Yr = HH*0;
for it=1:nrec
  tic;
  yr = YRPLT(it,1);
  mo = YRPLT(it,2);
  mday = YRPLT(it,3);
  dnmb = datenum(yr,mo,mday);
  iday = dnmb-datenum(yr,1,1)+1;
  fprintf('Extracting %i/%i/%i\n',yr,mo,mday);
  
  pthbin  = sprintf('/nexsan/archive/ARCc0.08_%3.3i/data/%4.4i_cice/',...
		    expt,yr);
  
  fin = sprintf('%s%3.3i_cice_inst.%4.4i-%2.2i-%2.2i-00000.nc',...
     	         pthbin,expt,yr,mo,mday);
%  keyboard
  hi = squeeze(nc_varget(fin,'hi')); % grid cell mean ice thickness
  voli = hi.*Acell;  % m3
%  ui = squeeze(nc_varget(fin,'uvel')); 
  vi = squeeze(nc_varget(fin,'vvel')); 
  
  
%  sflx = squeeze(nc_varget(fin,'fsalt_ai')); % S flx to ocean, kg/(m2 s)
%  fwflx= squeeze(nc_varget(fin,'fresh_ai'))*1e-2; %FW flux to ocean, m/mo
  
%  r=squeeze(nc_varget(fin,'rain_ai'))*1e-2; %rain, m/day
%  c=squeeze(nc_varget(fin,'congel'))*1e-2; %congelation, m/day
%  f=squeeze(nc_varget(fin,'frazil'))*1e-2; %frazil ice growth, m/day
%  mt=squeeze(nc_varget(fin,'meltt'))*1e-2;  %top ice melt, m/day
%  mb=squeeze(nc_varget(fin,'meltb'))*1e-2;  %bottom ice melt, m/day
%  ml=squeeze(nc_varget(fin,'meltl'))*1e-2;  %ltrl ice melt, m/day
%  s=squeeze(nc_varget(fin,'snow_ai'))*1e-2; %snowice form
%  fwf=r-c-f+mt+mb+ml+s;
%  fwf=r+mt+mb+ml+s;
%  pcolor(fwflx-fwf); shading flat;
%  caxis([-0.01 0.01]);
  
  cc=cc+1;
  for kb=1:nbx
% Estimate sea ice flux across bndry:
%    IJ=FIBX(kb).IJ;
%    i0=mean(IJ(:,1));
%    j0=mean(IJ(:,2));
    
    
    INP = FIBX(kb).IN_polygon;
    FIBX(kb).TM(cc,1)= dnmb;
    FIBX(kb).Title   = sprintf('Ice->ocean+ FW Flx m/mo, CICE-ARCc0.08-%i, %i/%i',expt,yr,mo);
 %   FIBX(kb).fwflx_m_mo(cc) = nansum(fwflx(INP).*Acell(INP))./sum(Acell(INP));
    FIBX(kb).ice_thkn(cc) = nansum(hi(INP).*Acell(INP))./sum(Acell(INP));
    FIBX(kb).ice_vol(cc) = nansum(hi(INP).*Acell(INP));
  end
  
% Calculate Fram Strait Ice Volume flux m3/s
  ifr1 = 930;
  ifr2 = 1078;
  jfr  = 950;

  vf=vi(jfr,ifr1:ifr2);
  hf=hi(jfr,ifr1:ifr2);
  dx=DX(jfr,ifr1:ifr2);
  fiflx = nansum(vf.*dx.*hf)*1e-9; % ice flux, km3/s
  FRIce(cc,1) = fiflx;

  
  fprintf('FW Sea Ice calculation %6.2f min\n',toc/60);
%  fprintf('%i: Beaufort Gyre, FW flux m/mo=%6.1f, Thkn m=%6.1fm\n',...
%	  yr,FIBX(7).fwflx_m_mo(cc),FIBX(7).ice_thkn(cc)); 
  fprintf('%i: Beaufort Gyre, FW Thkn m=%6.1fm, FramIceVolFlux=%6.1f km3/mo\n',...
	  yr,FIBX(7).ice_thkn(cc),fiflx*3600*24*30); 

  if s_mat>0 & mod(it,36)==0
    fprintf('Saving %s\n\n',fmat);
    save(fmat,'FIBX','FRIce');
  end;
end

% saving data
if s_mat>0
  fprintf('Saving %s\n\n',fmat);
  save(fmat,'FIBX','FRIce');
end;


