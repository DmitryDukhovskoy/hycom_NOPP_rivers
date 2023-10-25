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

mmean = 0; % = 1 - daily mean fields, 
regn = 'ARCc0.04';
%expt = 011;  
%expt = 012;  
expt = 022;

day_plot = datenum(2017,6,16);


% Test simulations with corrected CICE-HYCOM coupling
% Sept 2016, note different lenght of the runs
%expt = 0; % restart 0.04 GLBc GOFS3.5 created for 0.04 ARCc
%expt = 221; % fixed CICE-HYCOM coupled, restart from old run 0.04ARCc-CICEv5 Sept 1, 2016
%expt = 222; % fixed CICE-HYCOM coupled, CICE restart 0.04GLBc 01/09/2017
%expt = 223; % fixed CICE-HYCOM coupled, CICE&HYCOM restarts 0.04GLBc 01/09/2017 - CICE error T cond 9 days

plr   = 1;  % layer to plot
%plr   = 24;  % layer to plot, ~100m - similar to Claudia
pfld  = 'temp';
%pfld  = 'salin';

switch(pfld),
 case('salin');
%  c1=25.5;
%  c2=35.5;
   c1=29.5;
   c2=35.5;
%  c1=31;
%  c2=35.;
  ps=[1062 290 780 816];
 case('temp');
  c1=-2;
  c2= 1;
%  c2=12;
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
%for iyr=2017:2017
%  for idd=152:152  % 275 Oct 1, 336- Dec 1
%    if idd==1, idd=2; end;
%    cc=cc+1;
%    YRPLT(cc,1)=iyr;
%    YRPLT(cc,2)=idd;
%  end
%end

DV = datevec(day_plot);
iday = day_plot-datenum(2017,1,1)+1;
YRPLT(1,1) = DV(1);
YRPLT(1,2) = iday;

fprintf('Dates: %4.4i/%3.3i - %4.4i/%3.3i\n',YRPLT(1,1),YRPLT(1,2),...
	YRPLT(end,1),YRPLT(end,2));
np=size(YRPLT,1);


% Experiments:
% 102 - test simulation, HYCOM GLBb0.08 nest: 1 file kept constant (1993,1,1)

%pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/frames_TS/',regn,expt);
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
IDM=nn;
JDM=mm;

figure(1); clf;
set(gcf,'Position',[1209 280 1019 1055]);


% Plot fields:
cnc=0;
ip1=1;
for ip=ip1:np
  yr=YRPLT(ip,1);
  iday=YRPLT(ip,2);
%  pthbin = '/Net/mars/ddmitry/hycom/ARCc0.04/output/';
  switch(expt)
   case(12)
    pthbin = sprintf('/nexsan/hycom/ARCc0.04_%3.3i/data/%i/',expt,yr);  
    fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
    finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);
   case(22);
    pthbin = sprintf('/nexsan/people/ddmitry/hycom/ARCc0.04_022/data/%i/',yr);
    fina = sprintf('%s%3.3i_archv.%4.4i_%3.3i_00.a',pthbin,expt,yr,iday);
    finb = sprintf('%s%3.3i_archv.%4.4i_%3.3i_00.b',pthbin,expt,yr,iday);

    if mmean==1
      pthbin = sprintf('/nexsan/people/ddmitry/hycom/ARCc0.04_022/data/%i_mean/',yr);
      fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
      finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);
    end
   case(0); % restart
    pthbin = '/Net/kronos/ddmitry/hycom/ARCc0.04/restart_022/';
    fina = sprintf('%srestart_116i.a',pthbin);
    finb = sprintf('%srestart_116i.b',pthbin);
   case(221)
% Only last day of the run is saved
    yr=2016;
    iday=263;
    pthbin = '/nexsan/people/ddmitry/hycom/ARCc0.04_022/data_coupledCICE_restartOLDrun/';
    fina = sprintf('%s022_archv.%4.4i_%3.3i_00.a',pthbin,yr,iday);
    finb = sprintf('%s022_archv.%4.4i_%3.3i_00.b',pthbin,yr,iday);
   case(222)
    pthbin = '/Net/kronos/ddmitry/hycom/ARCc0.04/output_rest_ciceGOFS3.5_0901/';
    fina = sprintf('%sarchv.%4.4i_%3.3i_00.a',pthbin,yr,iday);
    finb = sprintf('%sarchv.%4.4i_%3.3i_00.b',pthbin,yr,iday);
   case(223)
    pthbin = '/Net/kronos/ddmitry/hycom/ARCc0.04/output_rest_allGOFS3.5_0901/';
    fina = sprintf('%s022_archv.%4.4i_%3.3i_00.a',pthbin,yr,iday);
    finb = sprintf('%s022_archv.%4.4i_%3.3i_00.b',pthbin,yr,iday);

  end
%  pthbin = sprintf('/nexsan/hycom/ARCc0.04_011/data012/%i/',yr);  % Greenland on exp
  
%[ZM,ZZ] = sub_zz_zm(fina,finb,HH);


  
  if ~exist(fina,'file');
    fprintf('Not found: %s\n\n',fina);
    continue;
  end
  
  cnc=cnc+1;
  dnmb=datenum(yr,1,1)+iday-1;
  DV=datevec(dnmb);


  fprintf('%s, sfig=%i, V.Lr=%i, %4.4i_%2.2i_%2.2i: %s\n',...
	  pfld,s_fig,plr,DV(1:3),fina);
% in archm u_vel=utot; total velocity
%  keyboard
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

 % Layer thickness:
  hfld = 'thknss';
  if expt==0, 
    hfld='dp'; 
    F = read_hycom_restart(fina,finb,hfld,IDM,JDM,'r_layer',plr); 
  else
    [F,n,m,l] = read_hycom(fina,finb,hfld,'r_layer',plr);
  end
  F=squeeze(F);
  F=F./rg;
  F(F>1e10)=0;
  F(F<1e-2)=0;
  dH=squeeze(F); 

  tic;
  if expt == 0   
    F = read_hycom_restart(fina,finb,pfld,IDM,JDM,'r_layer',plr);
  else
    [F,n,m,l] = read_hycom(fina,finb,pfld,'r_layer',plr);
  end

  F(F>1e6)=nan;

%  F = squeeze(F);
%  lTr = log(F);
  lTr = squeeze(F);

  nf = 1;
  ifx=max(strfind(fina,'/'));
  zLdp = mean_ZM_41lrs(plr);
  lTr(HH>zLdp) = nan;
  stl=sprintf('%s, %s, Lr %i, %6.1fm %i/%2.2i/%2.2i',...
	      fina(ifx+1:end),pfld,plr,zLdp,DV(1:3));
  
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
%  xlim1 = 1;
%  xlim2 = nn-1;
%  ylim1 = 1;
%  ylim2 = mm-1;

% Arctic Ocean:
  xlim1= 100;
  xlim2= nn;
  ylim1= 500;
  ylim2= 4000;



  if strncmp(pfld,'salin',4),
    f_cmp=6;
  else
    f_cmp=3;
  end

  clf;
  hps = [0.9 0.1 0.02 0.8];
  sub_plot_scalar(lTr,nf,HH,xlim1,xlim2,...
		  ylim1,ylim2,LON,LAT,stl,pfld,...
		  'c1',c1,'c2',c2,'cmp',f_cmp,'clbpos',hps);

  txtb='plot_arc04_scalar2.m';
  bottom_text(txtb,'pwd',1,'Position',[0.1 0.01 0.4 0.1]);

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





