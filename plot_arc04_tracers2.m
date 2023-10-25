% river passive tracers 
% distributed along the Greenland Coast
% Not depth-integrated 
% see anld_BG for depth-ingegrated plotting
% Tracer concentration in HYCOM experiment
% is set equal to actual FW flux (except for 
% Bering Strait) at the source = m3/sec
% In fact, this is max concentration that 
% the tracer can get after a relaxation time scale
% Thus, let's assume that the 
% simulated tracers = tracer concentration kg/m3
%
% Do not need to convert, as simulated is already concentration - kg/m3
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
s_reg = 0; %=0 - archm files are whole domain
         %=1 - archm files subsampled into nAtl. Region,
	 %      check fortran/hycom/extract_subdomain.F90
nTr=1;   % tracer to plot
TrMn=1e-10; % Threshold value to plot
rg = 9806;
fprintf('Tracer #: %i, Threhold value: %8.5d\n',nTr,TrMn);


% ======== BERING STRAIT ============
% Adjust Bering Strait tracer concenctration:
% Distribute FW flux (~0.08 Sv = 2500 km3/yr) over all points where Tracer=1
% # of j points = 14
% # of i points ~ 36
% total # of points (where HH>0 and above bottom)
% 438 points x 41 layers = 17958 grid points
% V.Flx = 79274 m3/s /17958 ~=4.4
% Thus tracer=1 unit is 4.4 m3/s of FW flux
% ---- This does not matter, assign some concentration
%      value, for comparison with other sources
%      set cBr = 15; 
cBr = 15.; % coefficient for Bering Str. tracer

f_zoom = 0;  % zoom in over the region of interest

%figure(1); clf;
%set(gcf,'Visible','off');
%figure('Visible','off'); clf;
%fprintf('Figure window is turned off\n');


%YPLT=[2004,299;2005,351;2006,295;2007,354;2008,286;2009,236];
% Animation:
YRPLT=[];
cc=0;
for iyr=2005:2005
  for idd=215:215
%for iyr=2008:2008
%  for idd=1:7:365
%    if idd==1, idd=2; end;
    cc=cc+1;
    YRPLT(cc,1)=iyr;
    YRPLT(cc,2)=idd;
  end
end

np=size(YRPLT,1);

% Experiments:
% 102 - test simulation, HYCOM GLBb0.08 nest: 1 file kept constant (1993,1,1)

regn = 'ARCc0.04';
expt = 011;  
%pthfig  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/fig_trac/',expt);
pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';
%pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc0.04/%3.3i/data_mat/',expt);

if s_reg==0 % full domain
  inc1=1;
  inc2=3200;
  jnc1=1;
  jnc2=5040;
  djnc=5040;
  dinc=3200;
elseif s_reg==1
%IND = smaller_domain_indices('Green');
  IND = subdomain_indx_arc04('atlarc');
  inc1=IND.i1;
  inc2=IND.i2;
  jnc1=IND.j1;
  jnc2=IND.j2;
  djnc=IND.dj;
  dinc=IND.di;
elseif s_reg==2
  inc1=250;  
  inc2=1590;
  jnc1=100;
  jnc2=2300;
end

ftopo = sprintf('%s/depth_%s_17DD.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');

HH=HH(jnc1:jnc2,inc1:inc2);
LON=LON(jnc1:jnc2,inc1:inc2);
LAT=LAT(jnc1:jnc2,inc1:inc2);
[mm,nn]=size(LON);

%hmsk=HH;
%hmsk(HH<0)=nan;

[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2

% Plot resolution:
f_prs=0;
if f_prs>0
  nf=5;
  xlim1=1; 
  xlim2=nn;
  ylim1=1;
  ylim2=mm;
  stl='Grid spacing, ARCc0.04';
  f_cmp=1;
  c1=1.75;
  c2=3.25;
  sub_plot_scalar_v3(DX*1e-3,nf,HH,xlim1,xlim2,...
		     ylim1,ylim2,LON,LAT,stl,c1,c2,f_cmp);
end;

xlim1 = 0;
xlim2 = nn;
ylim1 = 0;
ylim2 = mm;
if f_zoom
  switch(nTr)
   case(1)
    xlim1 = 100;
    xlim2 = 2200;
    ylim1 = 0;
    ylim2 = 2400;
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
%  pthbin = sprintf('/nexsan/archive/%s_%3.3i/data/%i/',regn,expt,yr)
%  pthbin = '/Net/mars/ddmitry/hycom/ARCc0.04/output/';
  pthbin = sprintf('/nexsan/hycom/ARCc0.04_011/data/%i/',yr);  

  if s_reg==0
    fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
    finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);
%    fina = sprintf('%s%3.3i_archv.%4.4i_%3.3i_00.a',pthbin,expt,yr,iday);
%    finb = sprintf('%s%3.3i_archv.%4.4i_%3.3i_00.b',pthbin,expt,yr,iday);
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

 % Layer thickness:
  [F,n,m,l] = read_hycom(fina,finb,'thknss','r_layer',plr);
  F=squeeze(F(plr,:,:));
  F=F./rg;
  F(F>1e10)=0;
  F(F<1e-2)=0;
  dH=squeeze(F); 

  tic;
  [F,n,m,l] = read_hycom(fina,finb,'tracer','r_tracer',nTr,'r_layer',plr);
  toc;
  if s_reg==0
    F=F(:,jnc1:jnc2,inc1:inc2);
  end
  F(F>1e6)=nan;
% For Pacific Water,
% Convert nondim concentration 1 to 
% actual FW flux-based conc.
% used Woodgate & Aagaard estimate of Bering FW flux
%
  if nTr==5, 
    fprintf('Adjusting Bering Strait flux\n');
    F = F*cBr;
  end
  Tr=squeeze(F(plr,:,:));

% Threshold value:
  Tr(Tr<=TrMn)=nan;
  I=find(~isnan(Tr));
% Do not need to convert, as simulated is already concentration - kg/m3
% Convert m3 of tr -> m3*1000 kg/m3/m3 x 1000(scaling) = kg(tracer)/m3(water)  
%  dmm = Tr(I)*1000./(dH(I).*Acell(I))*4000; % for plotting, scale x4000
%  Tr(I) = dmm;
  lTr=log(Tr);

  nf = 1;
%  stl=sprintf('Log2(C), Tr=%i, %4.4i/%2.2i/%2.2i, Layer %i',...
%	      nTr,DV(1:3),plr);
  ifx=max(strfind(fina,'/'));
  stl=sprintf('%s Log2(C, kg/m3), Tr=%i, Lr=%i',fina(ifx+1:end),nTr,plr);
  sub_plot_tracers(lTr,nf,HH,xlim1,xlim2,ylim1,ylim2,LON,LAT,stl,nTr);
  txtb='plot_arc04_tracers2.m';
  bottom_text(txtb,'pwd',1);
  
%keyboard
  
  %colorbar
%  drawnow
  
%  ffg=sprintf('%stracers_conc_%4.4i_%2.2i_%2.2i-layer%2.2i',pthfig,DV(1:3),plr);
%  ffg=sprintf('%s060trcr_lr%2.2i-%4.4i',pthfig,plr,cnc);
%  ffg=sprintf('%s060trcr_lr%2.2i-%3.3i',pthfig,plr,iday);

%  keyboard
  
  if s_fig>0
    txtb='plot_arc04_tracers2.m';
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




