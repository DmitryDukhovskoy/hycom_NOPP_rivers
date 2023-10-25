% Plot T or S fields from ARCc0.04 
% Monthly fields are extracted in extr_month_TS.m
%
% Plot only 1 layer at a time
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

s_fig = 0;

regn = 'ARCc0.04';
%expt = 011;  
expt = 022;  

plr   = 20;  % layer to plot
%plr   = 16;  % layer to plot, ~100m - similar to Claudia
pfld  = 'temp';
%pfld  = 'salin';

switch(pfld),
 case('salin');
%  c1=25.5;
%  c2=35.5;
  c1=31;
  c2=35.;
  ps=[1062 290 780 816];
 case('temp');
  c1=-2;
  c2=12;
  ps=[100 290 780 816];
end;  


pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';
%pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc0.04/%3.3i/data_mat/',expt);
pthout = '/Net/kronos/ddmitry/hycom/ARCc0.04/datamat/t_mnth/';
pthmat = pthout;


ftopo = sprintf('%s/depth_%s_17DD.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);


% Plot fields:
cnc=0;
ip1=1;
for yr=2016:2016
  for imo=8:8
    fmatout = sprintf('%s%s_monthly_lr%3.3i_%i%2.2i.mat',pthmat,pfld,plr,YR,im);

    if ~exist(fmatout,'file');
      fprintf('Not found: %s\n\n',fina);
      continue;
    end
  
% Greenland limits:
%  xlim1 = 1530;
%  xlim2 = 1900;
%  ylim1 = 1100;
%  ylim2 = 1600;
% Subpolar N. Atlantic
%  xlim1 = 600;
%  xlim2 = 2400;
%  ylim1 = 400;
%  ylim2 = 2200;
% Larbador Sea:
%  xlim1 = 650;
%  xlim2 = 1300;
%  ylim1 = 350;
%  ylim2 = 1400;
% Whole region
  xlim1 = 1;
  xlim2 = nn-1;
  ylim1 = 1;
  ylim2 = mm-1;

  if strncmp(pfld,'salin',4),
    f_cmp=6;
  else
    f_cmp=3;
  end

  clf;
  sub_plot_scalar(lTr,nf,HH,xlim1,xlim2,...
		  ylim1,ylim2,LON,LAT,stl,pfld,...
		  'c1',c1,'c2',c2,'cmp',f_cmp);

  txtb='plot_arc04_scalar.m';
  bottom_text(txtb,'pwd',1);

%  set(gcf,'position',ps);
  drawnow
  
  if s_fig>0
%    fnmF=sprintf('%s_lr%2.2i_%3.3i',pfld,plr,cnc);
    fnmF=sprintf('Labr_%s_lr%2.2i_%3.3i',pfld,plr,cnc);
    ffg=sprintf('%s%s',pthfig,fnmF);
    fprintf('Saving %s\n\n',ffg);
    print('-djpeg','-r200',ffg);
%    sso=sprintf('./trim_jpeg.com %s %s',pthfig,fnmF);
%    system(sso);
  end
  fprintf('Plotting 1 field: %6.4f min\n\n',toc/60);
  
end;  % day loop





