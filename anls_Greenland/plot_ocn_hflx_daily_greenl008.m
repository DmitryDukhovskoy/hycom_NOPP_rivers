% Plot ocean heat flux to Greenland
% across specified contour -
% isobath around Greenland
% saved daily flux - see 
% ocn_hflx_greenl008.m
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/seawater;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

YR1 = 2005;
YR2 = 2005;
Cp = 4200; % J/kg K
%Tref= -273.15; % Ref T to calc. H flux
Tref= -1.8; % Ref T to calc. H flux

expt = 110;


rg=9806;  % convert pressure to depth, m
hgg=1e20; 

plr=0; % highlight this interface
btx = 'plot_ocn_hflx_greenl008.m';

fprintf('Oceanic Heat Flux, Greenland Section, %i-%i\n',YR1,YR2);

regn = 'ARCc0.08';
expt = 110; % experiment without runoff
%expt = 112;  % epxeriment with Greenland runoff and monthly Arctic rivers
pthfig=sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/%s/%3.3i/fig_green_xsct/',...
		  regn,expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
%pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
pthmat  = '/Net/tholia/ddmitry/hycom/ARCc0.08/data_mat2/';

ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);

GC = sub_greenl_isobath(HH,LON,LAT);
Hs = GC.Hbottom; % bottom along section

% ================================================
% Plot Greenland map and the contour
% ================================================
f_pltgr=0;
if f_pltgr==1
  fn=10;
  sub_plot_Greenl_contour(HH,LON,LAT,fn,GC);
  
  bottom_text(btx,'pwd',1);
end


iyr=YR1;
fmat = sprintf('%s%3.3i_Greenl_HVflx_daily%i.mat',...
	       pthmat,expt,iyr);
load(fmat);

% plot total vol flux:
ndd=length(HFLX);
for i=1:ndd
  dmm=HFLX(i).Vol_flux_m3s;
end




  
  