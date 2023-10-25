% Plot mean T/S 
% GOFS3.1 reanalysis
% Expreiments: 53.X
% gridded netcdf 3hrly 
% extracted: extr_meanTS_GOFS31_GLB
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

s_fig  = 0;
pfld = 'salin';
zplt=0; 
%pfld = 'temp';
%mavg = 'DJF'; % winter/summer mean T/S
mavg = 'JJA'; % winter/summer mean T/S
yr1=2005;
yr2=2005;


rg=9806;  % convert pressure to depth, m
TV = '07';

pth190='/nexsan/GLBb0.08/GLBb0.08_190/data/meanstd/'; % 1993-1994
pth191='/nexsan/GLBb0.08/GLBb0.08_191/data/meanstd/'; % 1995-2012
pthmat = '/Net/mars/ddmitry/hycom/GLBb0.08/data_mat/';
pthfig = '/Net/mars/ddmitry/hycom/GLBb0.08/fig_fw/';
ftopo = sprintf('%sGLBb_T07_subset_Natl.mat',pthmat);
pthglb = '/nexsan/GLBb0.08/GLBb0.08_191/topo/';
%fmat = sprintf('%s%s_natl_%s.mat',pthmat,pfld,mavg);

fmat = sprintf('%sGOFS31_%s_natl_%4.4im_%i-%i_%s.mat',pthmat,pfld,abs(zplt),yr1,yr2,mavg);

txtb = 'plot_TS_GOFS31_GLB.m';

fprintf('Loading %s\n',fmat);
load(fmat);

AA = SAV.S_mean;
YR = SAV.YR;

%load(ftopo); doesn't match

switch(pfld);
  case('salin');
%   cl1=colormap_cold(200);
   cl1=flipud(colormap_cold(200));
   cmp=cl1;
   c1=31;
   c2=36;
 case('temp');
  nint = 200;
  c1=-2;
  c2=8;
  CMP = colormap_sclr2(nint,c1,c2);
  cmp = CMP.colormap;
  for k=1:15
    cmp(k,:)=[1 1 1];
  end
  cmp=smooth_colormap(cmp,15);
  cmp=smooth_colormap(cmp,15);

end

c1=30;
c2=35;

nint=200;
cmp = colormap(parula(nint));
cnt = (c1:(c2-c1)/nint:c2);


figure(1); clf;
pcolor(AA); shading flat;
hold on
axis('equal');
if strncmp(pfld,'salin',4);
  contour(AA,[30:0.5:36],'k');
  contour(AA,[33 33],'k','linewidth',1.6);
end

%xl1=300;
%xl2=930;
%yl1=10;
%yl2=800;
xl1=190;
xl2=962;
yl1=1;
yl2=800;
set(gca,'xlim',[xl1 xl2],...
	'ylim',[yl1 yl2],...
	'xticklabel',[],...
	'yticklabel',[],...
	'color',[0. 0. 0.]);

caxis([c1 c2]);
colormap(cmp);
hb=colorbar;
hps = [0.9 0.1 0.035 0.8];
set(hb,'position',hps,...
       'TickLength',0.04,...
       'Fontsize',22);

stl=sprintf('Global HYCOM+NCODA Reanlys. %s %s, %i-%i',pfld,mavg,YR);
title(stl);
btx = 'plot_TS_global_v2.m';
bottom_text(btx,'pwd',1);
