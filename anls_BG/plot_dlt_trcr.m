% Plot dlt concentration
% of tracers (relative to some year)
% for specifed year
% integrated over some depth:
% 0 - 50 m - mixed layer
% 150 - 50 m
%
% river passive tracers 
% distributed along the Greenland Coast
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

yr0 = 2005; 

lr1 = [0,-50];
lr2 = [-50,-150];


plr   = 1;  % layer to plot
s_fig = 1;
s_reg = 0; %=0 - archm files are whole domain
         %=1 - archm files subsampled into nAtl. Region,
	 %      check fortran/hycom/extract_subdomain.F90
nTr=5;   % tracer to plot
fprintf('Tracer #: %i\n',nTr);

f_zoom = 1;  % zoom in over the region of interest

%figure(1); clf;
%set(gcf,'Visible','off');
figure('Visible','off'); clf;
fprintf('Figure window is turned off\n');


%YPLT=[2004,299;2005,351;2006,295;2007,354;2008,286;2009,236];
% Animation:
YRPLT=[];
cc=0;
for iyr=2015:2015
  for idd=90:90
%for iyr=2008:2008
%  for idd=1:7:365
    if idd==1, idd=2; end;
    cc=cc+1;
    YRPLT(cc,1)=iyr;
    YRPLT(cc,2)=idd;
  end
end

np=size(YRPLT,1);

% Experiments:
% 102 - test simulation, HYCOM GLBb0.08 nest: 1 file kept constant (1993,1,1)

regn = 'ARCc0.08';
expt = 110;  
%pthfig  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/fig_trac/',expt);
pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc0.08/%3.3i/data_mat/',expt);

if s_reg==0
%IND = smaller_domain_indices('Green');
  IND = subdomain_indx_arc08('atlarc');
  inc1=IND.i1;
  inc2=IND.i2;
  jnc1=IND.j1;
  jnc2=IND.j2;
  djnc=IND.dj;
  dinc=IND.di;
elseif s_reg==1
  inc1=250;  
  inc2=1590;
  jnc1=100;
  jnc2=2300;
end

ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');

HH=HH(jnc1:jnc2,inc1:inc2);
LON=LON(jnc1:jnc2,inc1:inc2);
LAT=LAT(jnc1:jnc2,inc1:inc2);
[mm,nn]=size(LON);

hmsk=HH;
hmsk(HH<0)=nan;

xlim1 = 0;
xlim2 = nn;
ylim1 = 0;
ylim2 = mm;
if f_zoom
  switch(nTr)
   case(1)
    xlim1 = 50;
    xlim2 = 1100;
    ylim1 = 0;
    ylim2 = 1200;
  end
end




% Plot fields:
cnc=0;
ip1=1;
for ip=ip1:np
  yr=YRPLT(ip,1);
  iday=YRPLT(ip,2);
%  pthbin = sprintf(...
%    '/nexsan/people/ddmitry/hycom/ARCc0.08/%s/bin_outp/',expt);
  pthbin = sprintf('/nexsan/archive/%s_%3.3i/data/%i/',regn,expt,yr)

  if s_reg==0
    fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
    finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);
  elseif s_reg==1
    fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12nAtl.a',pthbin,expt,yr,iday);
    finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12nAtl.b',pthbin,expt,yr,iday);
  end
  
  if ~exist(fina,'file');
    fprintf('Not found: %s\n\n',fina);
    continue;
  end
  
  cnc=cnc+1;
  dnmb=datenum(yr,1,1)+iday-1;
  DV=datevec(dnmb);

%fprintf('Plotting: Tracer #%i V.Layer=%i\n\n',nTr,plr);
  fprintf('Tr# %i, sfig=%i, V.Lr=%i, %4.4i_%2.2i_%2.2i: %s\n',...
	  nTr,s_fig,plr,DV(1:3),fina);
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

  tic;
  [F,n,m,l] = read_hycom(fina,finb,'tracer','r_tracer',nTr,'r_layer',plr);
  toc;
  if s_reg==0
    F=F(:,jnc1:jnc2,inc1:inc2);
  end
  F(F>1e6)=nan;

%keyboard
  % Colormap:
  na=50;
  cl1=colormap_gray(na);
  cl2=colormap_purple(na);
  cl3=colormap_blue(na);
  cl4=colormap_green(na);
  cl5=colormap_yellow(na);
  cl6=colormap_red(na);

  cmp=[cl2;cl3;cl4;cl5;cl6];
  cmp=smooth_colormap(cmp,10);
  nint=length(cmp);
  c1=-6;
  c2=2;
  cnt=(c1:(c2-c1)/nint:c2);  % construct array of intervals


  if l>1
    Tr=squeeze(F(plr,:,:));
  else
    Tr=squeeze(F(plr,:,:));
  end
  
  Tr(Tr<=1e-3)=nan;
  lTr=log(Tr);
  % Normalize:
  %mT=max(max(Tr));
  %Tr=Tr./mT*100;
  clf
  pcolor(lTr); shading flat;
  hold on;
  contour(HH,[0 0],'k','Linewidth',1);
  contour(HH,[-5000:1000:-100],'Color',[0.7 0.7 0.7],'Linewidth',1);
  %caxis([0 2]);
  caxis([c1 c2]);
  colormap(cmp);
  axis('equal');
  set(gca,'xlim',[xlim1 xlim2],'ylim',[ylim1 ylim2]);
  set(gca,'xtick',[],'ytick',[]);
  clr=[0.9 0.9 0.9];
  plot_gridlines(45,10,1,clr,LON,LAT);
  stl=sprintf('Log2(C), Tr=%i, %4.4i/%2.2i/%2.2i, Layer %i',...
	      nTr,DV(1:3),plr);
  title(stl,'Fontsize',12);

  hght=[];
  lngth=[];
  mint=20;
  mbx=mint;
  fsz=12;
  bxc='k';
  posc=[0.81 0.1 0.8 0.06];
  aend=1;
  [az,axc]  = colorbar_vert(cmp,cnt,hght,lngth,mint,fsz,bxc,posc,mbx,aend);

  freezeColors;
  pcolor(hmsk); shading flat;
  colormap([0 0 0]);
  
  %colorbar
%  drawnow
  
%  ffg=sprintf('%stracers_conc_%4.4i_%2.2i_%2.2i-layer%2.2i',pthfig,DV(1:3),plr);
%  ffg=sprintf('%s060trcr_lr%2.2i-%4.4i',pthfig,plr,cnc);
%  ffg=sprintf('%s060trcr_lr%2.2i-%3.3i',pthfig,plr,iday);

%  keyboard
  
  if s_fig>0
    txtb='plot_tracers2.m';
    bottom_text(txtb,'pwd',1);
%    fnmF=sprintf('trcr%2.2i_lr%2.2i_%i-%4.4i',nTr,plr,yr,cnc);
    fnmF=sprintf('trcr%2.2i_lr%2.2i_%i-%3.3i',nTr,plr,yr,iday);
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




