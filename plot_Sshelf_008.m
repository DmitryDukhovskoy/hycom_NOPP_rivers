% Plot low S to highlight rivers
% sea ice edge September
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

plr   = 1;  % layer to plot
f_icec = 1; % plot saved monthly mean sea ice conc, NOAA
s_fig = 0;

% 110 - run with no Greenland, clim rivers
% 112 - run with Greenland runoff, NCAR rivers
expt = 112;  
%expt = 110;  
regn = 'ARCc0.08';
Fld  = 'salin';   % tracer to plot
fprintf('%s-%i, %s\n',regn,expt,Fld);


%YPLT=[2004,299;2005,351;2006,295;2007,354;2008,286;2009,236];
% Animation:
YRPLT=[];
cc=0;
for iyr=2010
  for idd=242:242
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
pthfig  = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_2d_flds/',regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc0.08/%3.3i/data_mat2/',expt);
%pthdat  = '/nexsan/people/ddmitry/Net_data2/ice_conc_GODDAR/monthly/';


ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

xlim1 = 270;
xlim2 = nn;
ylim1 = 580;
ylim2 = 1950;


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
  sub_plot_lowS(Fld,Tr,nf,HH,xlim1,xlim2,ylim1,ylim2,LON,LAT,stl);

% Ice conc extracted in
% get_iceconc_contours.m
  if f_icec==1
    fmatic=sprintf('%snoaa_iceconc_09_1987-1995.mat',pthmat);
    fprintf('Loading ice conc: %s\n',fmatic);
    load(fmatic);
    IC1=ICEC.Ice_edge_i;
    JC1=ICEC.Ice_edge_j;
    latc=ICEC.LAT;
    lonc=ICEC.LON;
    
    clear IH1 JH1
    npi=length(IC1);
    fprintf('Interpolating contour 1 to HYCOM grid ...\n');
    for ll=1:npi
      i0=IC1(ll);
      j0=JC1(ll);
      i1=floor(i0);
      j1=floor(j0);
      i2=i1+1;
      j2=j1+1;
      X=[i1 i2; i1 i2];
      Y=[j1 j1; j2 j2];
      V=[latc(j1,i1), latc(j1,i2); latc(j2,i1), latc(j2,i2)];
      lt0=interp2(X,Y,V,i0,j0);
      R=[lonc(j1,i1), lonc(j1,i2); lonc(j2,i1), lonc(j2,i2)];
      ln0=interp2(X,Y,R,i0,j0);
      
      d=dsphc(LAT,LON,lt0,ln0);
      [jm,im]=find(d==min(min(d)));
      IH1(ll,1)=im;
      JH1(ll,1)=jm;
    end
    
    fmatic=sprintf('%snoaa_iceconc_09_2010-2017.mat',pthmat);
    fprintf('Loading ice conc: %s\n',fmatic);
    load(fmatic);
    IC2=ICEC.Ice_edge_i;
    JC2=ICEC.Ice_edge_j;
    
    
    clear IH2 JH2
    npi=length(IC2);
    fprintf('Interpolating contour 2 to HYCOM grid ...\n');
    for ll=1:npi
      i0=IC2(ll);
      j0=JC2(ll);
      i1=floor(i0);
      j1=floor(j0);
      i2=i1+1;
      j2=j1+1;
      X=[i1 i2; i1 i2];
      Y=[j1 j1; j2 j2];
      V=[latc(j1,i1), latc(j1,i2); latc(j2,i1), latc(j2,i2)];
      lt0=interp2(X,Y,V,i0,j0);
      R=[lonc(j1,i1), lonc(j1,i2); lonc(j2,i1), lonc(j2,i2)];
      ln0=interp2(X,Y,R,i0,j0);
      
      d=dsphc(LAT,LON,lt0,ln0);
      [jm,im]=find(d==min(min(d)));
      IH2(ll,1)=im;
      JH2(ll,1)=jm;
    end
    
    plot(IH1,JH1,'r.');
    plot(IH2,JH2,'.','Color',[0 0.3 0.1]);

    axes('Position',[0.05 0.05 0.4 0.02]);
    sl{1}='Ice edge 0.15 conc, september 1987-1995(red), 2010-2017 (green)';
    sl{2}='Monthly sea ice conc, ';
    sl{3}='NOAA/NSIDC Climate Data Record of Passive Microwave';
    sl{4}='Sea Ice Concentration, Version 2 and 3';
    text(0,0,sl);
    set(gca,'visible','off');
  
  end

  
  txtb='plot_Sshelf_008.m';
  bottom_text(txtb,'pwd',1);
  
end;  % day loop





