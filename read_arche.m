% Read surface fluxes exchange btw ocean/sea ice
% output arche files:
%sst      =   sst, deg. C
%sss      =   sss
%ssu      =   ocean current, m/s
%ssv      =   
%ssh      =   sea surf. height, m
%ssfi     =   ocean heat flux to sea ice (dwnwrd), W/m2
%mlt      =   ocean mixed layer thickn, m
%sic      =   sea ice conc
%sitxdown =   downard sea ice tau x to ocean
%sitydown =   
%siqs     =   solar heat flux through ice to ocean, W/m2
%sifh     =   ice freezing/melting heat flux (-1 from CICE), W/m2
%sifs     =   ice freez/melt salt flux (kg / m2 s)
%sifw     =   ice net water flux, dwnrd (kg /m2 s)
%sit      =   sea ice T, (deg C)
%sih      =   sea ice thickness, m
%siu      =   sea ice x vel., m/s
%siv      =   
%surtx    =  surface wind flux (with no ice) 
%surty    =   
%sflice   =   ice salt flux

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

expt='022';
yr=2016;
iday=335;
hr=23;

j0=3200;  %Global i and j:        1354        1210
i0=5040;


%pthbin  = '/Net/kronos/ddmitry/hycom/ARCc0.08/tmp_files/';
%pthbin  = '/nexsan/people/ddmitry/hycom/ARCc0.08/110/surf/';
pthbin = '/nexsan/people/ddmitry/hycom/ARCc0.04_022/incoming/';
%pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';


%ftopo = sprintf('%sdepth_ARCc0.08_09.nc',pthtopo); % 
ftopo = sprintf('%s/depth_ARCc0.04_17DD.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');

%fina = sprintf('%s%s_arche.%4.4i_%3.3i_%2.2i.a',pthbin,expt,yr,iday,hr);
%finb = sprintf('%s%s_arche.%4.4i_%3.3i_%2.2i.b',pthbin,expt,yr,iday,hr);
fina = sprintf('%s022_arche.%4.4i_%3.3i_%2.2i.a',pthbin,yr,iday,hr);
finb = sprintf('%s022_arche.%4.4i_%3.3i_%2.2i.b',pthbin,yr,iday,hr);


%fld='ssfi';
fld='sitxdown'; % ice ocean stress
fld='surtx';    % wind ocean stress
[F,n,m,l] = read_hycom(fina,finb,fld);
F(F>1e20)=nan;
AA=squeeze(F);


f_plt = 1;
if f_plt==1
  figure(1); clf;
  pcolor(AA); shading flat;
  colorbar

%  caxis([-100 100]);
 
  idx = max(strfind(fina,'/'));
  stl = sprintf('%s, %s', fina(idx+1:end),fld);
  title(stl,'Interpreter','none');

end


