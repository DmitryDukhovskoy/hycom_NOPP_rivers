% Plot volume integrated mass of the tracer
% for specified regions
% extracted in vol_intgr_regn_trcr008.m
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

regn = 'ARCc0.08';
expt = 110;  
YR1 = 1993;
YR2 = 2016;

nTr   = 2;   % tracer to plot

% Experiments:
pthfig  = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc0.08/%3.3i/data_mat/',expt);
btx     = 'plot_VolIntgr_regn_trcr008.m';

%ftopo = sprintf('%s/depth_%s_17DD.nc',pthtopo,regn); % 
ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(HH);
[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2
MVOL = [];
Hmsk = HH;
%Hmsk(HH>h0)=100;


% Define Regions:
f_pltbox=0;  % plot boxes 
%BX = sub_define_boxes(HH,LON,LAT,f_pltbox);
%nbx = 5; % only N. Atl
BX = sub_define_boxesAO_Natl(HH,LON,LAT,f_pltbox);
for kk=1:length(BX)
  BX(kk).IN = BX(kk).IN_polygon;
end
%nbx = length(BX);
nbx=1; % Arctic Ocean

[XX,YY] = meshgrid((1:nn),(1:mm));

if ~isfield('BX','IN') 
  for ib=1:nbx
    iBG = BX(ib).IJ;
    INP = inpolygon(XX,YY,iBG(:,1),iBG(:,2));
    IN = find(INP==1);
    BX(ib).IN = IN;
  end
end



MVOL = struct;
ZM = [];
TM = [];
for ibx = 1:nbx
%  MVOL(ibx).mass_1m = [];
  MVOL(ibx).conc = [];
  MVOL(ibx).zm   = [];
end

for iyr=YR1:YR2
  yr   = iyr;
  fmat = sprintf('%strcr_regn_VrtLrMass_%i.mat',pthmat,iyr);
  
  fprintf('Loading %s\n',fmat);
  load(fmat);
  
%  tm = TRI.TM;
%  tm = tm(:)';
%  TM = [TM,tm];
% Mass tracers by layers:
  for ibx=1:nbx
    IN   = BX(ibx).IN;
    dzm  = VTRCR(nTr).DZM;
    zzm  = -cumsum(dzm,1);
    Mv   = VTRCR(nTr).TR; % kg
% Layers have different thickness, need to normalize to compare
% thus use concentration: kg/m3 for layer depths
    ar   = nansum(Acell(IN));
    Mvm  = Mv./(abs(dzm)*ar); % kg/m3
    ZM   = MVOL(ibx).zm;
    ZM   = [ZM,zzm];
%    mvol = MVOL(ibx).mass_1m;
    mvol = MVOL(ibx).conc;
    mvol = [mvol,Mvm];
    MVOL(ibx).conc = mvol; % kg/m3
    MVOL(ibx).zm      = ZM;
  end
  
%  keyboard
end; % time loop

cc=0;
for yr=YR1:YR2
  for im=1:12
    cc=cc+1;
    TM(cc)=datenum(yr,im,15);
  end
end


%DV  = datevec(TM);
%yrs = [DV(1,1):1/12:DV(end,1)+0.99];
yrs =[YR1:1/12:YR2+0.99];
[YRS,dmm]=meshgrid(yrs,ZM(:,1));

%nint=200;
%c1=0;
%c2=6;
%CMP = create_colormap5(nint,c1,c2);
%CMP = create_colormap4(nint,c1,c2);
%cnt = CMP.intervals;
%cmp = CMP.colormap;

nint = 320;
c1 = -8;
c2 = -5;
%CMP = create_colormap2_3(nint,c1,c2);
CMP = colormap_sclr2(nint,c1,c2);
cmp = CMP.colormap;
for ik=1:20
  cmp(ik,:) = [1 1 1];
end
cmp = smooth_colormap(cmp,20);
cmp = smooth_colormap(cmp,20);
cmp(1,:) = [1 1 1];
cnt = CMP.intervals;

% Plot annual vertical profiles for individual years 
%A=B
ibx = 1;
%itt = 288;
itt = 180;
Tr  = MVOL(ibx).conc(:,itt);
zm  = MVOL(ibx).zm(:,itt);
nm  = BX(ibx).Name;
dv= datevec(TM(itt));

figure(10+ibx); clf;
axes('Position',[0.2 0.15 0.4 0.8]);
plot(Tr,zm);
set(gca,'tickdir','out',...
	'ylim',[-3200 0],...
	'ytick',[-5000:500:0],...
	'ygrid','on',...
	'xgrid','on');

stl = sprintf('Tracer %i conc, %s, %4i/%2.2i',nTr,nm,dv(1:2));
title(stl);

bottom_text(btx,'pwd',1);

% Plot area-mean concentration
% i.e. total mass/1m divide by the box area - for comparison
for ibx = 1:nbx
% Box area:
  IN = BX(ibx).IN;
  amm = Acell(IN);
  Ar = nansum(amm);
  
  nm = BX(ibx).Name;
%  mvl = MVOL(ibx).mass_1m./Ar; % kg/1m -> kg/m3
  mvl = MVOL(ibx).conc; % kg/1m -> kg/m3
  ZM  = MVOL(ibx).zm;
  yl1 = 1.02*min(min(ZM));
  mvl(mvl<=0)=nan;
  lM  = log(mvl);

  figure(ibx); clf;
  axes('Position',[0.1 0.5 0.84 0.43]);
  pcolor(YRS,ZM,lM); shading interp;
  caxis([c1 c2]);
  set(gca,'tickdir','out',...
	  'xtick',[1990:2016],...
	  'xlim',[1993 2016],...
	  'ylim',[yl1 0],...
	  'ytick',[-3500:250:0]);
  
  colormap(cmp);
  
  
  stl = sprintf('%s, Tr %i, area-integrated',nm,nTr);
  title(stl);

  hc=colorbar('Location','SouthOutside');
  set(hc,'position',[0.1 0.42 0.84 0.012],...
      'Ticklength',0.014,'Fontsize',12);

  bottom_text(btx,'pwd',1,'position',[0.08 0.35 0.8 0.1]);
end




