% Spatial map of dltS for specified day
%
% Estimate S change in the upper layer due to surplus FW flux
% from Greenland 
% Note estimate of Greenland runoff fraction
% for specified boxes may use wrong
% code - tracer concentration averaged over layers
% should be integrated giving a total mass within
% a layer
% the Fraction is calculate in fraction_tracer_NAtl.m
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

LR(1,:) = [0,-50];
LR(2,:) = [-50,-150];
LR(3,:) = [0,-10000]; % full water column
LR(4,:) = [-150, -300];
nlrs= 4; 

regn = 'ARCc0.08';
expt = 110;  
pthfig  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/fig_NAtl/',expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
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


pthbin  = sprintf('/nexsan/archive/ARCc0.08_110/data/%4.4i/',yr);
fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);
fld='salin';
[SS,n,m,l] = read_hycom(fina,finb,fld);
SS(SS>1e20)=nan;
[ZM,ZZ] = sub_zz_zm(fina, finb,HH,'f_btm',1);



BX = sub_define_boxes(HH,LON,LAT,0);

[XX,YY] = meshgrid((1:nn),(1:mm));
for ib=1:5
  iBG = BX(ib).IJ;
  INP = inpolygon(XX,YY,iBG(:,1),iBG(:,2));
  IN = find(INP==1);
  BX(ib).IN = IN;
end

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
hZ=abs(HH);
hZ(HH>=0)=nan;
ilv=3;
rrD = sub_fraction_tracerMass(expt,regn,Acell,HH,dnmb,nTr,ilv,hZ);


for ilv=1:4
  zz1 = LR(ilv,1);
  zz2 = LR(ilv,2);

  dz=abs(zz2-zz1);

  fprintf('S change in regions , ilv=%i, %i - %i, nTr=%i, %i/%i/%i\n',ilv,zz1,zz2,nTr,dv0(1:3));


  % Find mean S for layer ilv
  Zk  = HH*0;
  cZk = HH*0; % count v. layers
  for kk=1:l
    zz= squeeze(ZZ(kk,:,:));
    I=find(HH<zz2 & zz>=zz2 & zz<zz1);
    if ilv==3
      zzn= squeeze(ZZ(kk+1,:,:));
      I=find(~isnan(zz) & ~isnan(zzn)); % find bottom
    end
    if isempty(I); continue; end;
    if max(max(zz))<zz2, break; end;
    Zk(I)=kk;
    cZk(I)=cZk(I)+1;
  end

  Zk(Zk==0)=nan; % bottom layer within the depth interval
  Zt=Zk-cZk+1;   % top layer
  Sav = HH*0; % S averaged within the layer
  lmx = max(max(Zk));
  hZ  = HH*0;
  smm = HH*0;
  for kk=1:lmx
    I=find(Zk>=kk & Zt<=kk);
    if isempty(I), continue; end;
    dZ=abs(squeeze(ZZ(kk+1,:,:))-squeeze(ZZ(kk,:,:)));
    hZ(I)=hZ(I)+dZ(I);
  %  if isnan(hZ(300,600)); keyboard; end;
    smm = smm+squeeze(SS(kk,:,:)).*dZ;
    Sav(I) = smm(I)./hZ(I); % depth-avrg S
  end
  Sav(Sav==0)=nan;

% Tracer fraction in grid cells
  rr = sub_fraction_tracerMass(expt,regn,Acell,HH,dnmb,nTr,ilv,hZ);

% Estimate volume of Greenland surplus FW in grid cell
% = ratio of FW tracer*total FWF (integrated over time 1990-date requested)
  Vfw = fwf0*rr*1e9; % m3 
  Hfw = Vfw./Acell; % m of FW

% Estiamte dlt(S) due to Greenland FWFlux anomaly
  Vol = Acell.*hZ;
  dS=-(Sav-(Vol-Vfw).*Sav./Vol); % freshening -> should be negative
  dS(abs(dS)<1e-23)=nan;

  
% By regions
  for ib=1:5
    nm=BX(ib).Name;
    if ilv==1, xlb{ib}=nm; end;
    
    IN=BX(ib).IN;
    dmm=nanmean(Hfw(IN));
    FWC(ib,ilv)=dmm;
    dmm=nanmean(dS(IN));
    DS(ib,ilv)=dmm;
  end
end

% Estimate dltS if all FW stays in the upper layers: 50, 150, 300 m
% Note here S change is for the whole layer from surface
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

  for ib=1:5
    nm=BX(ib).Name;
    IN=BX(ib).IN;
    dmm=nanmean(HfwD(IN));
    FWCD(ib,kk)=dmm;
    dmm=nanmean(dSD(IN));
    DSD(ib,kk)=dmm;
  end
end




dmm=FWC;
FWC(:,3)=dmm(:,4);
FWC(:,4)=dmm(:,3);
dmm=DS;
DS(:,3)=dmm(:,4);
DS(:,4)=dmm(:,3);

% Delete Whole depth data
FWC  = FWC(:,1:3);
DS   = DS(:,1:3);

% Bar Diagram:
figure(1); clf;
axes('Position',[0.08 0.55 0.88 0.4]);
bar(FWC,0.95);
sll=sprintf('FWC (m) change, 2016');
title(sll);
set(gca,'xticklabel',xlb,...
	'tickdir','out');
%lg=legend('0-50','50-150','150-300','btm');
lg=legend('0-50','50-150','150-300');

axes('Position',[0.08 0.08 0.88 0.4]);
bar(DS,0.95);
sll=sprintf('dS, 2016');
title(sll);
set(gca,'xticklabel',xlb,...
	'tickdir','out');

btx='dSregions_bars.m';
bottom_text(btx,'pwd',1)


cmp = colormap_blue(nl);
% FWC and dS if all Gr. FW in the region
% stays in the layer
% This shows dS for the whole layer surface to dZ
figure(2); clf;
axes('Position',[0.08 0.55 0.88 0.4]);
hb = bar(DSD,0.95);
for ik=1:nl
  clr=cmp(ik,:);
  hb(ik).FaceColor=clr;
end

sll=sprintf('dS FW stays in layers, 2016');
title(sll);
set(gca,'xticklabel',xlb,...
	'tickdir','out');
colormap(cmp);
hc=colorbar('horizontal');
set(hc,'TickLength',0.018);

btx='dSregions_bars.m';
bottom_text(btx,'pwd',1,'position',[0.08 0.5 0.6 0.1])



