% Extract time series of monthly 
% mean area vorticity
% in the BG to estimate strength of BG
% similar to AOO
% divergence is calculated in ../anls_UV/calc_divU_month.m
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers
startup;

close all
clear

yr1=1993;
yr2=2016;
zz0 =50;

regn = 'ARCc0.08';
expt = 110;  
%pthfig  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/fig_trac/',expt);
pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = '/Net/tholia/ddmitry/hycom/ARCc0.08/data_mat2/';
btx = 'calc_mean_vort_ocean.m';

ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);
Acell = DX.*DY;

BX  = sub_define_boxesAO_Natl(HH,LON,LAT,0);
nbg = 7; 
IN  = BX(nbg).IN_polygon;

cc=0;
for iyr = yr1:yr2
%  fmat = sprintf('%s%3.3i_divU%3.3im_%i.mat',...
%		   pthmat,expt,abs(zz0),iyr);
  fmat = sprintf('%s%3.3i_divU%3.3im_v2_%i.mat',...
		   pthmat,expt,abs(zz0),iyr);
  fprintf('Loading %s\n',fmat);
  load(fmat);
  
  for imo=1:12
    cc=cc+1;
    dvu = DIVU(imo).divU_m3_sec;
    Vrt(cc,1) = nanmean(dvu(IN)./(Acell(IN)*abs(zz0))); % area mean vort, s-1
    fprintf('%i/%2.2i mean vort: %6.4d\n',iyr,imo,Vrt(cc));
  end
end

yrs=[yr1:1/12:yr2+0.99];
for yr=yr1:yr2
  I=find(yrs>=yr & yrs<yr+1);
  Vanmn(I)=mean(Vrt(I));
end

figure(1); clf;
axes('Position',[0.08 0.45 0.85 0.45]);
plot(yrs,Vrt);
hold on
plot([yr1 yr2+1],[0 0],'-','Color',[0.6 0.6 0.6]);
plot(yrs,Vanmn,'r-','Linewidth',1.5);
set(gca,'tickdir','out',...
	'xlim',[yr1 yr2+1],...
	'ylim',[-3.5e-8 3.5e-8],...
	'xtick',[yr1:yr2],...
	'ytick',[-3e-8:1e-8:3e-8],...
	'yminortick','on',...
	'xminortick','on',...
	'xgrid','on',...
	'ygrid','on',...
	'Fontsize',14);

title('Area-Mean Ocean U vorticity, BG, s-1');

set(gcf,'Position',[622 479 1870 863]);

bottom_text(btx,'pwd',1,'Position',[0.05 0.3 0.8 0.1],'Fontsize',9);






  
