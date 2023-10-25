% Plot 2D 
% distributed along the Greenland Coast
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

plr   = 1;  % layer to plot
s_fig = 0;

% 110 - run with no Greenland, clim rivers
% 112 - run with Greenland runoff, NCAR rivers
%expt = 112;  
expt = 110;  
regn = 'ARCc0.08';
Fld  = 'salin';   % tracer to plot
fprintf('%s-%i, %s\n',regn,expt,Fld);


%YPLT=[2004,299;2005,351;2006,295;2007,354;2008,286;2009,236];
% Animation:
YRPLT=[];
cc=0;
for iyr=1995
  for idd=230:230
%for iyr=2008:2008
%  for idd=1:7:365
    if idd==1, idd=2; end;
    cc=cc+1;
    YRPLT(cc,1)=iyr;
    YRPLT(cc,2)=idd;
  end
end
np=size(YRPLT,1);

%pthfig  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/fig_trac/',expt);
pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_2d_flds/',regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc0.08/%3.3i/data_mat/',expt);


ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

xlim1 = 0;
xlim2 = nn;
ylim1 = 0;
ylim2 = mm;

% Plot fields:
cnc=0;
ip1=1;
for ip=ip1:np
  yr=YRPLT(ip,1);
  iday=YRPLT(ip,2);

  switch(expt),
   case(112)
    pthbin = sprintf('/nexsan/hycom/ARCc0.08_112/data/%i/',yr);
   case(110),
    pthbin = sprintf('/nexsan/archive/ARCc0.08_110/data/%4i/',yr);
  end

  fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
  finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);
  
  if ~exist(fina,'file');
    fprintf('Not found: %s\n\n',fina);
    continue;
  end
  
  cnc=cnc+1;
  dnmb=datenum(yr,1,1)+iday-1;
  DV=datevec(dnmb);

%fprintf('Plotting: Tracer #%i V.Layer=%i\n\n',nTr,plr);
  fprintf('%s, sfig=%i, V.Lr=%i, %4.4i_%2.2i_%2.2i: %s\n',...
	  Fld,s_fig,plr,DV(1:3),fina);
% in archm u_vel=utot; total velocity
  f_getu=0;
  if f_getu==1
    ilr=1;
    [Fu,n,m,l] = read_hycom(fina,finb,'u-vel.');
    [Fv,n,m,l] = read_hycom(fina,finb,'v-vel.');
    U=squeeze(Fu(ilr,:,:));
    U(U>1e6)=nan;
    V=squeeze(Fv(ilr,:,:));
    V(V>1e6)=nan;
    S=sqrt(U.^2+V.^2);
  end

  tic;
  [F,n,m,l] = read_hycom(fina,finb,Fld,'r_layer',plr);
  toc;
  F(F>1e6)=nan;

  if l>1
    Tr=squeeze(F(plr,:,:));
  else
    Tr=squeeze(F(plr,:,:));
  end
  
  Tr(Tr<=1e-3)=nan;

%keyboard
  nf = 1;
  ifx=max(strfind(fina,'/'));
  stl=sprintf('%s, %s, Lr=%i, %i/%i/%i',fina(ifx+1:end),Fld,plr,DV(1:3));
  sub_plot_2dfld(Fld,Tr,nf,HH,xlim1,xlim2,ylim1,ylim2,LON,LAT,stl);
  txtb='plot_archm2D_008.m';
  bottom_text(txtb,'pwd',1);
  
%keyboard
  
  %colorbar
%  drawnow
  
%  ffg=sprintf('%stracers_conc_%4.4i_%2.2i_%2.2i-layer%2.2i',pthfig,DV(1:3),plr);
%  ffg=sprintf('%s060trcr_lr%2.2i-%4.4i',pthfig,plr,cnc);
%  ffg=sprintf('%s060trcr_lr%2.2i-%3.3i',pthfig,plr,iday);

%  keyboard
  
  if s_fig>0
%    txtb='plot_archm2D_008.m';
%    bottom_text(txtb,'pwd',1);
%    fnmF=sprintf('trcr%2.2i_lr%2.2i_%i-%4.4i',nTr,plr,yr,cnc);
    fnmF=sprintf('%s_lr%2.2i_%i-%3.3i',Fld,plr,yr,iday);
    ffg=sprintf('%s%s',pthfig,fnmF);
    fprintf('Saving %s\n\n',ffg);
    print('-djpeg','-r200',ffg);
%    sso=sprintf('./trim_jpeg.com %s %s',pthfig,fnmF);
%    system(sso);
  end
end;  % day loop

%I=find(~isnan(Tr));
%T=Tr(I);
%b=[1e-30;1e-10;1e-6;1e-3;1e-2;1e-1;1;(10:100)'];
%N=hist(T,b);
%bar(log10(b),N);




