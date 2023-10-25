% Calculate 
% Tracer fluxes across POP gates
% on Greenland shelf - extracted in greenl_Trcr_POPgates.m
%
% The gate coordinates are in the following order: Davis Strait, 
% Cape Farewell, Wide to Narrow, Denmark Strait, Fram Strait, and  Narres Strait. 

addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/Wavelet;
startup;

close all
clear


expt=110;
TV=11;
YR1=2013;
YR2=2016;


pthmat = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_theresa/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';

btx='plot_greenl_Trcr_POPgates.m';

ftopo = sprintf('%sdepth_ARCc0.08_%2.2i.nc',pthtopo,TV); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);

GC = sub_greenl_isobath(HH,LON,LAT);
Ig=GC.cntr_Iindx;
Jg=GC.cntr_Jindx;

GrFWA_dn=[];
GrFWA_dv=[];
jjy=0;
for YR=YR1:YR2
%  if YR==2009; continue; end
  jjy=jjy+1;
  yr=YR;
  dE=datenum(yr,12,31);
  dJ1=datenum(yr,1,1);
  
  fmatout=sprintf('%shycom008_%3.3i_Greenl_Trcr_POPgates_%4.4i.mat',...
		pthmat,expt,YR);
  fprintf('Loading %s\n',fmatout);
  load(fmatout);

  Ftr_dv=VTRCR(1).TrFlx;  % tracer flux, kg/s, over the whole section
  Fsh_dv=VTRCR(1).TrFlx_Sh;
  Vtr_dv=VTRCR(1).VolFlx_m3s;
  Vsh_dv=VTRCR(1).VolFlxGrSh_m3s;

  Ftr_dn=VTRCR(4).TrFlx;  % tracer flux, kg/s, over the whole section
  Fsh_dn=VTRCR(4).TrFlx_Sh;
  Vtr_dn=VTRCR(4).VolFlx_m3s;
  Vsh_dn=VTRCR(4).VolFlxGrSh_m3s;

  dday=7;
  TM=[dJ1:dday:dJ1+364]; 
  DVm=datevec(TM);

% Scale everything to tracer mass  
% First, get overall tracer mass in the whole domain:
  pthmat2=sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
  nTr=1;

  % Get Greenland FW anomaly for given date:
% Read in cumulative Greenland FW flux anomaly
% mat file created in dSregions_map_v2.m
  frv = sprintf('%sGreenland_cumFWFlux.mat',pthmat2);
  fprintf('Loading %s\n',frv);
  load(frv);

  imo0=-1;
  for it=1:length(TM)
    DV = DVm(it,:);
    iyr = DV(1);
    imo = DV(2);
    if imo~=imo0
      fmat2 = sprintf('%sMassTr%2.2i_lrs_%i%2.2i.mat',pthmat2,nTr,iyr,imo);
      fprintf('Loading %s\n',fmat2);
      load(fmat2);
    end
% find whole-depth layer
    nlr = length(TRCR);
  %if isfield('Layer_Z2',TRCR);
    ibtm=5; % whole-depth tracer mass
%  for ilr=1:nlr
%    z2=TRCR(ilr).Layer_Z2;
%    if z2<-9000, ibtm=ilr; end;
%  end
    Tr_dom=squeeze(TRCR(ibtm).MassTr_kg); % tracer integrated ove whole water depth
    amm=Tr_dom./(abs(HH).*DX.*DY);  % kg/m3
    Tr_dom(amm<=0.125)=nan;

  %vol = Acell*hZ; % grid cell vol
  %intgrTr_dom = nansum(nansum(Tr.*Acell.*hZ));
    MTr_dom = nansum(nansum(Tr_dom));

% Greenland runoff anomaly for given date
    ism = find(Ygr==DV(1,1));
    fwf0 = cFWF(ism); % km3

% Convert tracer mass flux kg/s to GFW anomaly flux km3/mo
% Average by months
    MTr_dom(MTr_dom==0)=nan;
    rr=Ftr_dv(it)/MTr_dom;
    rr2=Ftr_dn(it)/MTr_dom;
    vflx=rr*fwf0*3600*24*30; % km3/mo
    vflx2=rr2*fwf0*3600*24*30; 
    if imo~=imo0 | it==length(TM)
      if imo0>0
        GrFWA_dv(jjy,imo0)=dmm/icc; % GFWA anomaly flux cross Davis strait, km3/mo
        GrFWA_dn(jjy,imo0)=dmm2/icc; % Denmark Str, km3/mo
      end
      icc=0;
      imo0=imo;
      dmm=0;
      dmm2=0;
    end
    icc=icc+1;
    dmm=dmm+vflx;
    dmm2=dmm2+vflx2;
  end

end


%mF1=mean(GrFWA_dv);
%mF2=mean(GrFWA_dn);

mF1=median(GrFWA_dv);
p11=prctile(GrFWA_dv,90);
p12=prctile(GrFWA_dv,10);
mF2=median(GrFWA_dn);
p21=prctile(GrFWA_dn,90);
p22=prctile(GrFWA_dn,10);

clr1=[0 0.4 0.8];
clr2=[0.8 0.4 0];

figure(1); clf;
axes('Position',[0.1 0.6 0.85 0.32]);
hold on;
pp1=plot(mF1,'-','Linewidth',2,'Color',clr1);
plot(mF1,'.','Markersize',14,'Color',clr1);

pp2=plot(mF2,'-','Linewidth',2,'Color',clr2);
plot(mF2,'.','Markersize',14,'Color',clr2);

for ik=1:12
  y1=p11(ik);
  y2=p12(ik);
  plot([ik ik],[y1 y2],'--','color',clr1,'Linewidth',1);
  plot([ik-0.05 ik+0.05],[y1 y1],'-','color',clr1,'Linewidth',1);
  plot([ik-0.05 ik+0.05],[y2 y2],'-','color',clr1,'Linewidth',1);
 
  y1=p21(ik);
  y2=p22(ik);
  plot([ik ik],[y1 y2],'--','color',clr2,'Linewidth',1);
  plot([ik-0.05 ik+0.05],[y1 y1],'-','color',clr2,'Linewidth',1);
  plot([ik-0.05 ik+0.05],[y2 y2],'-','color',clr2,'Linewidth',1);
end

ym1=1.05*floor(min([min(min(mF1)),min(min(mF2)),min(p12),min(p22)]));
ym2=1.05*ceil(max([max(max(mF1)),max(max(mF2)),max(p11),max(p21)]));

set(gca,'xlim',[0.9 12.1],...
        'ylim',[ym1 ym2],...
        'xtick',[1:12],...
        'ytick',[-50:5:50],...
        'tickdir','out',...
        'Fontsize',14,...
        'xgrid','on',...
        'ygrid','on');

lgg=legend([pp1,pp2],'Davis','Denmark');
set(lgg,'Position',[0.6 0.35 0.24 0.11],'Fontsize',14);
xlabel('Months');
ylabel('GFWA, km3/mo');
ttl=sprintf('HYCOM 0.08-%3.3i, GFWA Transport, %i-%i',expt,YR1,YR2);
title(ttl);

bottom_text(btx,'pwd',1,'Position',[0.05 0.3 0.4 0.05]);

%
%
% Time integrated flux of GFWA -> km3/yr at the end

cF1=cumsum(mF1);
cF2=cumsum(mF2);

figure(2); clf;
axes('Position',[0.1 0.6 0.85 0.32]);
hold on;

plot(cF1,'-','Linewidth',2,'Color',clr1);
plot(cF1,'.','Markersize',14,'Color',clr1);
plot(cF2,'-','Linewidth',2,'Color',clr2);
plot(cF2,'.','Markersize',14,'Color',clr2);

ym1=1.05*floor(min([min(cF1),min(cF2)]));

ym2=1.05*ceil(max([max(cF1),max(cF2)]));
ym2=max([ym1, ym2, 0]);

set(gca,'xlim',[0.9 12.1],...
        'ylim',[ym1 ym2],...
        'xtick',[1:12],...
        'ytick',[-70:5:50],...
        'tickdir','out',...
        'Fontsize',14,...
        'xgrid','on',...
        'ygrid','on');

lgg=legend([pp1,pp2],'Davis','Denmark');
set(lgg,'Position',[0.6 0.35 0.24 0.11],'Fontsize',14);
xlabel('Months');
ylabel('GFWA, km3');
ttl=sprintf('HYCOM 0.08-%3.3i, Cumulative GFWA Transp, Median %i-%i',expt,YR1,YR2);
title(ttl);

bottom_text(btx,'pwd',1,'Position',[0.05 0.3 0.4 0.05]);



