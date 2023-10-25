% Plot mean T/S 
% Greenland - see extr_meanTS_Glb.m
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

s_fig  = 0;
pfld = 'salin';
%pfld = 'temp';
%mavg = 'DJF'; % winter/summer mean T/S
mavg = 'JJA'; % winter/summer mean T/S


rg=9806;  % convert pressure to depth, m
TV = '07';

pth190='/nexsan/GLBb0.08/GLBb0.08_190/data/meanstd/'; % 1993-1994
pth191='/nexsan/GLBb0.08/GLBb0.08_191/data/meanstd/'; % 1995-2012
pthmat = '/Net/mars/ddmitry/hycom/GLBb0.08/data_mat/';
pthfig = '/Net/mars/ddmitry/hycom/GLBb0.08/fig_fw/';
ftopo = sprintf('%sGLBb_T07_subset_Natl.mat',pthmat);
pthglb = '/nexsan/GLBb0.08/GLBb0.08_191/topo/';
fmat = sprintf('%s%s_natl_%s.mat',pthmat,pfld,mavg);
txtb = 'plot_TS_global.m';

fprintf('Loading %s\n',fmat);
load(fmat);

AA = SAV.S_mean;
YR = SAV.YR;

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

figure(1); clf;
pcolor(AA); shading flat;
hold on
axis('equal');
if strncmp(pfld,'salin',4);
  contour(AA,[30:0.5:36],'k');
  contour(AA,[33 33],'k','linewidth',1.6);
end

set(gca,'xlim',[300 830],...
	'ylim',[100 800],...
	'xticklabel',[],...
	'yticklabel',[],...
	'color',[0.4 0.4 0.4]);

caxis([c1 c2]);
colormap(cmp);
hb=colorbar;
set(hb,'Fontsize',22);

stl=sprintf('Global HYCOM+NCODA Reanlys. %s %s, %i-%i',pfld,mavg,YR);
title(stl);
btx = 'plot_TS_global.m';
bottom_text(btx,'pwd',1);
