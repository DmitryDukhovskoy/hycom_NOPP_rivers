% Plot winter (ONDJFM) mean heat fluxes OAFlux
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

s_fig  = 0;


rg=9806;  % convert pressure to depth, m
Tref=-1.9; % Reference T, C
Zref=-500; % depth for Heat Cont. calc. 
Cp = 4186; % heat capacity, J*kg/C
TV = '07';

pth190='/nexsan/GLBb0.08/GLBb0.08_190/data/meanstd/'; % rean. 1993-1994
pth191='/nexsan/GLBb0.08/GLBb0.08_191/data/meanstd/'; % rean. 1995-2012
pthmat = '/Net/mars/ddmitry/hycom/GLBb0.08/data_mat/';
pthdata = '/Net/yucatan/tachanat/flux_products/';
pthfig = '/Net/mars/ddmitry/hycom/GLBb0.08/fig_HeatCont/';
ftopo = sprintf('%sGLBb_T07_subset_Natl.mat',pthmat);
pthglb = '/nexsan/GLBb0.08/GLBb0.08_191/topo/';

txtb = 'ddmitry/.../hycom_NOPP_rivers/validation_GLBb/plot_AOFlux.m';

YRPLT=[];
cnc=0;
for yr=2000:2015
  cnc=cnc+1;
  YRPLT(cnc,1)=yr;
end

fprintf('Loading topo %s\n',ftopo);
load(ftopo);
[mm,nn]=size(HH);

% Grid cell spacing
[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2

HFLUX = sub_define_regions(HH);
NR=length(HFLUX);

fmatin=sprintf('%sCERES_OAFlux_Qnet_2000_2015_monthly.mat',pthdata);
F=load(fmatin);
lonf=F.lon;
I=find(lonf>180);
lonf(I)=lonf(I)-360;
latf=F.lat;
[X,Y]=meshgrid(lonf,latf);
nnf=length(lonf);
mmf=length(latf);
[DXf,DYf]=sub_dx_dy(X,Y);
AreaF=DXf.*DYf; % Grid cell area, m2


for ik=1:NR
  IJs=HFLUX(ik).IJs;
  i1=min(IJs(:,1));
  i2=max(IJs(:,1));
  j1=min(IJs(:,2));
  j2=max(IJs(:,2));
  ln1=LON(j1,i1);
  lt1=LAT(j1,i1);
  ln2=LON(j2,i2);
  lt2=LAT(j2,i2);
  dd=abs(lonf-ln1);
  ii1=find(dd==min(dd),1);
  dd=abs(lonf-ln2);
  ii2=find(dd==min(dd),1);
  dd=abs(latf-lt1);
  jj1=find(dd==min(dd),1);
  dd=abs(latf-lt2);
  jj2=find(dd==min(dd),1);
  
  MM=zeros(mmf,nnf);
  if ii2>=ii1
    MM(jj1:jj2,ii1:ii2)=1;
  else
    MM(jj1:jj2,ii1:nnf)=1;
    MM(jj1:jj2,1:ii2)=1;
  end
  IN=find(MM==1);
  HFLUX(ik).INflx=IN;
end

npl=length(YRPLT);
for ipp=2:npl
  yr1=YRPLT(ipp-1);
  yr2=YRPLT(ipp);
  
  for ir=1:NR
    HFLUX(ir).Qsum=0; % area averaged
    HFLUX(ir).Qall=[]; % percentiles
    HFLUX(ir).Winter_Year(ipp-1)=yr2;
  end
  
  cc=0;
  for im=12:12  % months
    cc=cc+1;
    fnml=sprintf('dmm=F.qnet%4im%i;',yr1,im);
    eval(fnml);
    for ir=1:2
      IN=HFLUX(ir).INflx;
      qnet=nansum(dmm(IN).*AreaF(IN))./nansum(AreaF(IN)); %area averaged
      Qall=HFLUX(ir).Qall;
      Qall=[Qall;dmm(IN)];
      HFLUX(ir).Qall=Qall;  % 90th percentile for all winter months
      HFLUX(ir).Qsum=HFLUX(ir).Qsum + qnet;
    end
  end

  for im=1:2 % months
    cc=cc+1;
    fnml=sprintf('dmm=F.qnet%4im%i;',yr2,im);
    eval(fnml);
    for ir=1:NR
      IN=HFLUX(ir).INflx;
      qnet=nansum(dmm(IN).*AreaF(IN))./nansum(AreaF(IN));
      Qall=HFLUX(ir).Qall;
      Qall=[Qall;dmm(IN)];
      HFLUX(ir).Qall=Qall;  % 90th percentile for all winter months
      HFLUX(ir).Qsum=HFLUX(ir).Qsum + qnet;
    end
  end
  for ir=1:NR
    HFLUX(ir).Qnet_wint(ipp-1,1)=HFLUX(ir).Qsum/cc;
    Qall=HFLUX(ir).Qall;
    q90=prctile(Qall,90);  % 90th percentile heat flux over the region
    HFLUX(ir).Qperc90(ipp-1,1)=q90;
  end
  
  
end

POS = [0.05 0.5 0.35 0.4;
       0.05 0.05 0.35 0.5];
CLR=[0,0,0; 0,0,1; 0,1,0; 1,0,0];

%figure(1); clf;
%axes('Position',[0.08 0.55 0.88 0.35]);
hold on;
for ir=1:NR
  qflx = HFLUX(ir).Qnet_wint;
  yrp  = HFLUX(ir).Winter_Year;
  clr=CLR(ir,:);
%  plot(yrp,qflx,'Linewidth',2,'Color',clr);
  Qnet(:,ir)=qflx;
end
title('Winter Mean Net Heat Flux, + to atm');

%axes('Position',[0.08 0.05 0.88 0.35]);
hold on;
for ir=1:NR
  qflx = HFLUX(ir).Qperc90;
  yrp  = HFLUX(ir).Winter_Year;
  clr=CLR(ir,:);
%  plot(yrp,qflx,'Linewidth',2,'Color',clr);
  Q90p(:,ir)=qflx;
  NMS{ir}=HFLUX(ir).Name;
end
title('90prc Winter Net Heat Flux, + to atm');

cmp=[0,0.5,.5; 0,1,1; 0,1,0; 0.8,0.8,0];

figure(2); clf;
axes('Position',[0.08 0.55 0.88 0.35]);
hb=bar(yrp,Q90p,0.95);
colormap(cmp);
hlg=legend(NMS);
set(hlg,'Position',[0.65 0.28 0.11 0.14]);
set(gca,'Tickdir','out',...
	'ylim',[0 310],...
	'xlim',[2000.5 2015.5],...
	'xtick',[2001:2015],...
	'Fontsize',14);

sttl='AOFlux Monthly Mean, 90prcntile of DJF';
title(sttl,'Fontsize',14);
bottom_text(txtb);

if s_fig>0
  fgnm=sprintf('%sAOFlx_winter_90prcnt_regions_2001-2015',pthfig);
  fprintf('Saving %s\n',fgnm);
  print('-dpng','-r250',fgnm);
end


%axes('Position',[0.08 0.05 0.88 0.35]);
f_reg=0;
if f_reg>0
% Plot regions
  nf=10;
  nrplt=[];
  sub_plot_regions(nf,HFLUX,HH,LON,LAT,nrplt);
  bottom_text(txtb);
  
  if s_fig>0
    fgnm=sprintf('%sAOFlx_regions',pthfig);
    fprintf('Saving %s\n',fgnm);
    print('-dpng','-r250',fgnm);
  end

end
    
