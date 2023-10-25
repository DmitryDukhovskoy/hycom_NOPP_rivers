% Bar diagrams dltS for specified day
% Updated code
% uses tracer integrated over specified v. layers
% integration uses exact depth levels
% extr_MassTrcr_month.m
% 
% Also uses depth- & month-averaged S as a reference S
% to calculate FW gain and dltS
%
% Estimate S change in the upper layer due to surplus FW flux
% from Greenland 
%
% Tracer content calculated in vol_intgr_regn_trcr008.m
%
% Note: updated regions, now all regions
% are adjacent to each other
% So that Subpolar Gyre
% combines reg #1 (Labr)+Reg#2(IcelSea)+Reg#7 (Centr NAtl)
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers
startup;

close all
clear

s_fig = 0;

dnmb = datenum(2016,12,30);
dv0  = datevec(dnmb);
yr   = dv0(1);
iday = dnmb-datenum(yr,1,1)+1;

nTr = 1; 
%ilv = 1; % 0-50m
%ilv = 2; % 50-150m
%ilv = 3; % bottom-surf
%ilv = 4; % 150-300m


% Vertical layers
LRS = load('LRS.dat');
nlrs= length(LRS); 

regn = 'ARCc0.08';
expt = 110;  
pthfig  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/fig_NAtl/',expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
pthsav  = sprintf('/Net/tholia/ddmitry/hycom/ARCc0.08/data_mat%3.3i/',expt);
%fmout   = sprintf('%sfraction_trcr_NAtl.mat',pthmat);

ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

% Grid cell spacing
[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2


hmsk=HH;
hmsk(HH<0)=nan;

% Old regions
%BX = sub_define_boxes(HH,LON,LAT,0);
% New regions SPNA
h0  = -800; % cutoff depth, m
f_pltbox=0;
BX = sub_deep_regions(HH,LON,LAT,h0,f_pltbox);
%BX = sub_define_boxes(HH,LON,LAT,f_pltbox);
%nbx = length(BX);
nbx = length(BX); % 

[XX,YY] = meshgrid((1:nn),(1:mm));
%for ib=1:5
%  iBG = BX(ib).IJ;
%  INP = inpolygon(XX,YY,iBG(:,1),iBG(:,2));
%  IN = find(INP==1);
%  BX(ib).IN = IN;
%end

f_griv = 0; %=1 - rederive cumulative Greenland FWFlux from Bamber
frv = sprintf('%sGreenland_cumFWFlux.mat',pthmat);
if f_griv==1
  [TMg,Fgr]=sub_read_Greenland_v3; % Greenland total FWF, km3/mo
  % Calculate annual anomalies
  % and cumulative up to date
  DVg = datevec(TMg);
  ngr = length(TMg);
  nyr =ngr/12;
  dmm = reshape(Fgr,[12,nyr]);

  Fyr = sum(dmm); % km3/yr
  Ygr = [DVg(1,1):DVg(end,1)];
  ii=find(Ygr==1990);
  Fmn = mean(Fyr(1:ii));
  ism = find(Ygr==dv0(1,1));
  cFWF = cumsum(Fyr-Fmn);
  fwf0 = cFWF(ism); % km3
  save(frv,'cFWF','Ygr');
else
  fprintf('f_griv %i, Loading %s\n',f_griv,frv);
  load(frv);
%  ii=find(Ygr==1990);
  ism = find(Ygr==dv0(1,1));
  fwf0 = cFWF(ism); % km3
end  


% Whole Depth:
% Tracer fraction in grid cells
%hZ=abs(HH);
%hZ(HH>=0)=nan;
%ilv=3;
%rrD = sub_fraction_tracerMass(expt,regn,Acell,HH,dnmb,nTr,ilv,hZ);


for ilv=1:4
  zz1 = LRS(ilv,1);
  zz2 = LRS(ilv,2);
  dz  = abs(zz2-zz1);

  fprintf('S change in regions , ilv=%i, %i - %i, nTr=%i, %i/%i/%i\n',ilv,zz1,zz2,nTr,dv0(1:3));


%
% Get layer-averaged S - as a reference S
% to calculate dlt S
% prepared in anls_TS/mnthly_arc08_layers_S.m
  YR=dv0(1);
  mo=dv0(2);
  fsout = sprintf('%sarc08_%3.3i_mnthS_lrs_%4.4i%2.2i.mat',pthsav,expt,YR,mo);
  fprintf('Loading S averaged %s\n',fsout);
  load(fsout);
  Sav = meanS(ilv).Savrg;
  Sav(Sav==0)=nan;
  hZ = abs(LRS(ilv,2)-LRS(ilv,1));
  
% Tracer fraction in grid cells
  IntSrf=0;
  rr = sub_fraction_tracerMass(expt,regn,Acell,HH,dnmb,nTr,ilv,hZ,IntSrf);

% Estimate volume of Greenland surplus FW in grid cell
% = ratio of FW tracer*total FWF (integrated over time 1990-date requested)
  Vfw = fwf0*rr*1e9; % m3 
  dZrf = 50; % refence water column height
  Hfw = Vfw./Acell*dZrf/dz; % m of GrFW in dZrf m of water column

% Estiamte dlt(S) due to Greenland FWFlux anomaly
  Vol = Acell.*hZ;
  dS=-(Sav-(Vol-Vfw).*Sav./Vol); % freshening -> should be negative
  dS(abs(dS)<1e-23)=nan;

  
% By regions
  for ib=1:nbx
    nm=BX(ib).Name;
    if ilv==1, xlb{ib}=nm(1:6); end;
    
    IN=BX(ib).IN;
    dmm=nanmean(Hfw(IN));
    FWC(ib,ilv)=dmm;
    dmm=nanmean(dS(IN));
    DS(ib,ilv)=dmm;
  end
end

% Estimate dltS if all FW stays in the upper layers: 50, 150, 300 m
% Note here S change is for the whole layer from surface
ilv = 5; % whole depth
IntSrf=0;
rrD = sub_fraction_tracerMass(expt,regn,Acell,HH,dnmb,nTr,ilv,hZ,IntSrf);
LRS=[50:50:500];
nl=length(LRS);
clear DSD FWCD
for kk=1:nl
  dz=LRS(kk);
  Vol=dz.*Acell;
  VfwD = fwf0*rrD*1e9; % m3, FW vol integrated over the whole Depth
  HfwD = VfwD./Acell; % m of FW
  dSD  = -Sav+(Vol-VfwD).*Sav./Vol; % freshening negative
  dSD(abs(dSD)<1e-23)=nan;

  for ib=1:nbx
    nm=BX(ib).Name;
    IN=BX(ib).IN;
    dmm=nanmean(HfwD(IN));
    FWCD(ib,kk)=dmm;
    dmm=nanmean(dSD(IN));
    DSD(ib,kk)=dmm;
  end
end


%dmm=FWC;
%FWC(:,3)=dmm(:,4);
%FWC(:,4)=dmm(:,3);
%dmm=DS;
%DS(:,3)=dmm(:,4);
%DS(:,4)=dmm(:,3);

% Delete Whole depth data
%FWC  = FWC(:,1:3);
%DS   = DS(:,1:3);

cmp = [0.98 0.98 0.98;...
       0.8 0.8 0.8;...
       0.5 0.5 0.5;...
       0.2 0.2 0.2];

% Bar Diagram:
figure(1); clf;
axes('Position',[0.08 0.55 0.88 0.4]);
hb=bar(FWC,0.95); 
for ik=1:4
  set(hb(ik),'FaceColor',cmp(ik,:));
end
sll=sprintf('GrFWC (m per %i m), 2016',dZrf);
title(sll);
set(gca,'xticklabel',xlb,...
	'tickdir','out');
%lg=legend('0-50','50-150','150-300','btm');
lg=legend('0-50','50-150','150-300','300-500');

axes('Position',[0.08 0.08 0.88 0.4]);
hb=bar(DS,0.95);
for ik=1:4
  set(hb(ik),'FaceColor',cmp(ik,:));
end
sll=sprintf('dS, 2016');
title(sll);
set(gca,'xticklabel',xlb,...
	'tickdir','out',...
	'Fontsize',14);

btx='dSregions_bars_v2.m';
bottom_text(btx,'pwd',1)


%cmp = colormap_blue(nl);
cmp = colormap_gray(nl);
% FWC and dS if all Gr. FW in the region
% stays in the layer
% This shows dS for the whole layer surface to dZ
figure(2); clf;
axes('Position',[0.08 0.43 0.88 0.5]);
hb = bar(DSD,0.95);
for ik=1:nl
  clr=cmp(ik,:);
  hb(ik).FaceColor=clr;
end
yl1 = 1.1*min(min(DSD));

sll=sprintf('dS if GrFW stays in layers, 2016');
title(sll);
set(gca,'xticklabel',xlb,...
	'tickdir','out',...
	'ylim',[yl1 0],...
	'Fontsize',14);
colormap(cmp);
hc=colorbar('horizontal');
set(hc,'TickLength',0.018);

bottom_text(btx,'pwd',1,'position',[0.08 0.35 0.6 0.1])



