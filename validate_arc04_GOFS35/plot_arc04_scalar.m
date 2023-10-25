% Plot T or S fields from ARCc0.04 
% Plot only 1 layer at a time
% if need to layer-average - modify the code
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

rg = 9806;
if s_fig==1;
  fprintf('Plotting field, %s fig saved ON, Layer: %i\n',pfld,plr);
else
  fprintf('Plotting field, Not saved: %s, Layer: %i\n',pfld,plr);
end

YRPLT=[];
cc=0;
for iyr=2016:2016
  for idd=32:32
    if idd==1, idd=2; end;
    cc=cc+1;
    YRPLT(cc,1)=iyr;
    YRPLT(cc,2)=idd;
  end
end

fprintf('Dates: %4.4i/%3.3i - %4.4i/%3.3i\n',YRPLT(1,1),YRPLT(1,2),...
	YRPLT(end,1),YRPLT(end,2));
np=size(YRPLT,1);

% Experiments:
% 102 - test simulation, HYCOM GLBb0.08 nest: 1 file kept constant (1993,1,1)

pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';
%pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc0.04/%3.3i/data_mat/',expt);

%inc1=1;
%inc2=1600;
%jnc1=1;
%jnc2=2520;
%djnc=2520;
%dinc=1600;

ftopo = sprintf('%s/depth_%s_17DD.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2


% Plot fields:
cnc=0;
ip1=1;
for ip=ip1:np
  yr=YRPLT(ip,1);
  iday=YRPLT(ip,2);
%  pthbin = '/Net/mars/ddmitry/hycom/ARCc0.04/output/';
%  pthbin = sprintf('/nexsan/hycom/ARCc0.04_%3.3i/data/%i/',expt,yr);  
%  pthbin = '/Net/kronos/ddmitry/hycom/ARCc0.04/output/';
  pthbin = '/Net/kronos/ddmitry/hycom/ARCc0.04/output_eloan/';

%  pthbin = sprintf('/nexsan/hycom/ARCc0.04_011/data012/%i/',yr);  % Greenland on exp
  fina = sprintf('%sarchv.%4.4i_%3.3i_00.a',pthbin,yr,iday);
  finb = sprintf('%sarchv.%4.4i_%3.3i_00.b',pthbin,yr,iday);
  fina = sprintf('%s022_archv.%4.4i_%3.3i_00.a',pthbin,yr,iday);
  finb = sprintf('%s022_archv.%4.4i_%3.3i_00.b',pthbin,yr,iday);
  
  
  if ~exist(fina,'file');
    fprintf('Not found: %s\n\n',fina);
    continue;
  end
  
  cnc=cnc+1;
  dnmb=datenum(yr,1,1)+iday-1;
  DV=datevec(dnmb);


  fprintf('%s, sfig=%i, V.Lr=%i, %4.4i_%2.2i_%2.2i: %s\n',...
	          pfld,s_fig,plr,DV(1:3),fina);
 
% Layer thickness:
  [F,n,m,l] = read_hycom(fina,finb,'thknss','r_layer',plr);
  F=squeeze(F);
  F=F./rg;
  F(F>1e10)=0;
  F(F<1e-2)=0;
  dH=squeeze(F); 

  tic;
  [F,n,m,l] = read_hycom(fina,finb,pfld,'r_layer',plr);

  F(F>1e6)=nan;

%  F = squeeze(F);
%  lTr = log(F);
  lTr = squeeze(F);

  nf = 1;
  ifx=max(strfind(fina,'/'));
%  zLdp = mean_ZM_41lrs(plr);
%  lTr(HH>zLdp) = nan;
  stl=sprintf('%s, %s, Lr %i, %i/%2.2i/%2.2i',...
	      fina(ifx+1:end),pfld,plr,DV(1:3));
  
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





