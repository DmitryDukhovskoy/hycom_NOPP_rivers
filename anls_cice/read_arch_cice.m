% Read CICE output
%

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

expt='110';
yr=1993;
mo=1;
dm=2;

%pthbin  = '/Net/kronos/ddmitry/hycom/ARCc0.08/tmp_files/';
%pthbin  = '/nexsan/people/ddmitry/hycom/ARCc0.08/110/cice/';
pthbin = sprintf('/nexsan/archive/ARCc0.08_%s/data/%i_cice/',expt,yr);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';

ftopo = sprintf('%sdepth_ARCc0.08_11.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');

fin = sprintf('%s%s_cice_inst.%4.4i-%2.2i-%2.2i-00000.nc',...
	      pthbin,expt,yr,mo,dm);

%A=squeeze(nc_varget(fin,'divu'));
%u=squeeze(nc_varget(fin,'uocn'));
%v=squeeze(nc_varget(fin,'vocn'));
%u=squeeze(nc_varget(fin,'uvel'));
%v=squeeze(nc_varget(fin,'vvel'));
%S=sqrt(u.^2+v.^2);

%sairx=squeeze(nc_varget(fin,'strairx')); % N/m2
%sairy=squeeze(nc_varget(fin,'strairy'));
%S=sqrt(sairx.^2+sairy.^2);

%pcolor(S); shading flat

hi = squeeze(nc_varget(fin,'hi'));
pcolor(hi); shading flat;
caxis([0 3]);
colorbar




