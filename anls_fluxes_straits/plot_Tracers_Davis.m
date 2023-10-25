% Plot snapshot of Tracer Conc (Arctic only, no Greenland Tr.) 
% across Fram Strait
% 
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear


dnmb = datenum(2016,12,1);
dv = datevec(dnmb);
YR = dv(1);
mo = dv(2);
mday = dv(3);
yr = YR;
iday = dnmb-datenum(yr,1,1)+1;

expt = 110;
ntr  = 4;
hgg  = 1e20;
rg   = 9806;


pthmat  = '/Net/mars/ddmitry/hycom/ARCc0.08/data_mat/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
btx = 'plot_Tracers_Davis.m';

cBr = sub_BeringTrConc(0);% coefficient for Bering Str., tracer, kg/m3

fprintf('Tracer Conc, Davis Strait  %s\n',datestr(dnmb));

ftopo = sprintf('%sdepth_ARCc0.08_11.nc',pthtopo); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);

%SCT = sub_set_Fram_xsct(HH,LON,LAT);
SCT = sub_set_Davis_xsct(HH,LON,LAT);


pthbin = sprintf('/nexsan/archive/ARCc0.08_%3.3i/data/%4.4i/',expt,yr);
fina   = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
finb   = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);
  
   
    
tic;
%    [F,n,m,nlr] = read_hycom(fina,finb,'temp');
%    F(F>hgg)=nan;
%    T=F;

%    [F,n,m,nlr] = read_hycom(fina,finb,'salin');
%    F(F>hgg)=nan;
%    S=F;
%
%[F,n,m,nlr] = read_hycom(fina,finb,'u-vel.');
%F(F>hgg)=0;
%U=F;

[F,n,m,nlr] = read_hycom(fina,finb,'v-vel.');
F(F>hgg)=0;
V=F;

%fld='thknss';
%[F,n,m,l] = read_hycom(fina,finb,fld);
%F(F>1e18)=0;
%F=F/rg;
%F(F<1e-2)=0;
%dH = F;
[ZM,ZZ] = sub_zz_zm(fina, finb, HH);

% Process segments
II = SCT.II;
JJ = SCT.JJ;
np = length(II);

c1=-6;
c2=4;
nint=200;
cnt=(c1:(c2-c1)/nint:c2);  % construct array of intervals
Rr=-1;                        % colors of max intensity, clockwise rotation 
C0=[-1,0,1];                 % starting point, red dimension is Z axis
Cend=[1, 1, 0.6];
cmpT = colormap_spiral(nint,'C0',C0,'Rr',Rr,'Cend',Cend);
for k=1:15
  cmpT(k,:)=[1 1 1];
end
cmpT=smooth_colormap(cmpT,10);

% ----------------
% Tracrers
% ----------------
TRN{1}='Mackenzie';
TRN{2}='EEuras';
TRN{3}='WEuras';
TRN{4}='Bering';
for nTr = 2:5
  fprintf('   Reading Tr# %i\n',nTr);

  [F,nn,mm,ll] = read_hycom(fina,finb,'tracer','r_tracer',nTr);
  F(F>1e6)=nan;

% For Pacific Water,
% Convert nondim concentration 1 to 
% actual FW flux-based conc.
% used Woodgate & Aagaard estimate of Bering FW flux
  if nTr==5, 
    F = F*cBr;
  end
  F(F<=0)=nan;
  Ctr = log(F); % kg/m3
  
  j1=JJ(1);
  i1=II(1);
  i2=II(end);
  vv=squeeze(V(:,j1,i1:i2));
  ctr=squeeze(Ctr(:,j1,i1:i2));
  Hb=HH(j1,i1:i2);
  zz=squeeze(ZZ(:,j1,i1:i2));
  xx=LON(j1,i1:i2);
  [XX,dmm]=meshgrid(xx,zz(:,1));
  ym=mean(LAT(j1,i1:i2));
  
  ctr=[ctr;ctr(end,:)];
  vv=[vv;vv(end,:)];
  
  figure(nTr); clf;
  axes('Position',[0.09 0.1 0.8 0.8]);
  pcolor(xx,zz,ctr); shading interp
  hold on
  caxis([c1 c2]);
  colormap(cmpT);
  contour(XX,zz,vv,[0:0.05:0.8],'k-','linewidth',1.);
  contour(XX,zz,vv,[-0.8:0.05:-0.001],'k--','linewidth',1.);
  
  vxb=[xx(1),xx,xx(end)];
  vyb=[-10000,Hb,-10000];
  patch(vxb,vyb,[0 0 0]);
 
  stl=sprintf('008-%3.3i, Davis,%3.1fN, log(Tr) %s, %i/%i/%i',...
	      expt,ym,TRN{nTr-1},dv(1:3));
  title(stl,'Fontsize',14);
  
  set(gca,'tickdir','out',...
	  'xlim',[xx(1) xx(end)],...
	  'ylim',[-1000 0],...
	  'ytick',[-4000:200:0],...
	  'Fontsize',14);
  
  hb=colorbar;
  set(hb,'Position',[0.91 0.2 0.015 0.6],...
	 'TickLength',0.02,...
	 'Ticks',[c1:1:c2],...
	 'Fontsize',16);

  set(gcf,'Position',[1019 594 1289 748]);
  
  bottom_text(btx,'pwd',1);  
  
  
end
% ----------------
fprintf('1 day: %8.4f min\n\n',toc/60);
    
