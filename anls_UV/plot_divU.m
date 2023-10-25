% Plot mean div U calculated
% in calc_divU_month.m
%

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers
startup;

close all
clear

yr1=2000;
yr2=2016;

zz0 = 50;

regn = 'ARCc0.08';
expt = 110;  
%pthfig  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/fig_trac/',expt);
pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_mat2/';
btx = 'plot_divU.m';

ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);
Acell = DX.*DY;

% Filtering mask
pgrd = 51;
Hmsk = HH;
Hmsk(HH<0)=1;
Hmsk(HH>=0)=0;
Hmsk(1920:end,:)=0;
Hmsk(100:500,1250:end)=0;
Hmsk(1:pgrd+1,:) = 0;
Hmsk(mm-pgrd:mm,:) = 0;
Hmsk(:,1:pgrd+1) = 0;
Hmsk(:,nn-pgrd:nn) = 0;

cl2 = colormap_red(100);
cl1 = colormap_blue(100);
for ik=1:2;
  cl2(ik,:) = [1 1 1];
  cl1(ik,:) = [1 1 1];
end
cl1 = flipud(cl1);
cmp = [cl1;cl2];
cmp = smooth_colormap(cmp,10);

cc=0;
DVU=HH*0;
for iyr = yr1:yr2
%  fmat = sprintf('%s%3.3i_divU%3.3im_%i.mat',...
%		   pthmat,expt,abs(zz0),iyr);
  fmat = sprintf('%s%3.3i_divU%3.3im_v2_%i.mat',...
		   pthmat,expt,abs(zz0),iyr);
  fprintf('Loading %s\n',fmat);
  load(fmat);
  
  for imo=1:12 % bug, Dec is missing
    fprintf('mo=%i\n',imo);
    cc=cc+1;
    dvu = DIVU(imo).divU_m3_sec;
    DVU=DVU+dvu;
%    dvuf = sub_fltr(dvu,pgrd,Hmsk);
    
  end
  
  

end

DVU=DVU/cc;

figure(imo); clf;
pcolor(dvu); shading flat;
colormap(cmp);
%caxis([-5e4 5e4]);
caxis([-1e3 1e3]);

hold on;
contour(HH,[0 0],'k');

colorbar
title(sprintf('DivU, m3/s, upper 50m, %i-%i',yr1,yr2));
axis('equal');
set(gca,'xlim',[320 1100],...
	'ylim',[300 1100],...
	'color',[0 0 0]);

bottom_text(btx,'pwd',1);



  



