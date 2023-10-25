% Plot annual heat content in the N. Atl
% Plot HC anomaly relaitve to 1990s mean
% prepared in heat_cont_natl_annual.m
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

s_fig  = 1;


rg=9806;  % convert pressure to depth, m
Tref=-1.9; % Reference T, C
Zref=-100; % depth for Heat Cont. calc. 
Cp = 4186; % heat capacity, J*kg/C
TV = '07';

pth190='/nexsan/GLBb0.08/GLBb0.08_190/data/meanstd/'; % rean. 1993-1994
pth191='/nexsan/GLBb0.08/GLBb0.08_191/data/meanstd/'; % rean. 1995-2012
pthmat = '/Net/mars/ddmitry/hycom/GLBb0.08/data_mat/';
pthfig = '/Net/mars/ddmitry/hycom/GLBb0.08/fig_HeatCont/';
ftopo = sprintf('%sGLBb_T07_subset_Natl.mat',pthmat);
pthglb = '/nexsan/GLBb0.08/GLBb0.08_191/topo/';
%fmat = sprintf('%sHeatCont_natl_1993-2015.mat',pthmat);
fmat = sprintf('%sHeatCont%4.4im_natl_1993-2015.mat',pthmat,abs(round(Zref)));
txtb = 'ddmitry/.../hycom_NOPP_rivers/validation_GLBb/plot_annual_heatcont.m';

YRPLT=[];
cnc=0;
for yr=1993:2015
  cnc=cnc+1;
  YRPLT(cnc,1)=yr;
end

fprintf('Loading topo %s\n',ftopo);
load(ftopo);
[mm,nn]=size(HH);

% Grid cell spacing
[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2


fprintf('Loading heat content data %s\n',fmat);
load(fmat);

nyrs=length(HTC);
for iy=1:nyrs
  YR(iy,1)=HTC(iy).Year;
end

Zref=HTC(1).Zref;
Tref=HTC(1).Tref;
cff=1e-9; % J -> GJ

% Calculate mean:
htc0=zeros(mm,nn);
cc=0;
for iy=1:nyrs
  cc=cc+1;
  htc=HTC(iy).HeatCnt_J_m2*cff; % delta htc
%    htc(htc<1e-1)=nan;
  htc0=htc0+htc;
end
htc0=htc0/cc;



%YRPLT=[1993; 1998; 2008; 2015];
YRPLT=[];
npl=length(YRPLT);
if npl>0
  for iyr=1:npl
%    htc0=HTC(1).HeatCnt_J_m2*cff;
    year=YRPLT(iyr);
    iy=find(YR==year);
%    htc=HTC(iy).HeatCnt_J_m2*cff;
    htc=HTC(iy).HeatCnt_J_m2*cff-htc0; % delta htc
%    htc(htc<1e-1)=nan;

    figure(iyr); clf;
    hold on;
    contour(HH,[0 0],'k');
    pcolor(htc); shading flat;
    caxis([-4 4])
    colorbar
%    ttl=sprintf('GLBb0.08, HeatCont%6.0d J/m2, intigr. Z=%5.1f,Tref=%4.1f,  %i',...
%		1/cff,Zref,Tref,year);
    ttl=sprintf('GLBb0.08, dlt(HeatCont)%6.0d to 1993 J/m2, Z=%5.1f,Tref=%4.1f,  %i',...
		1/cff,Zref,Tref,year);
    title(ttl,'Fontsize',14);
    txtb='ddmitry/.../hycom_NOPP_rivers/validation_GLBb/plot_annual_heatcont.m';
    bottom_text(txtb);
    
    set(gca,'xtick',[],...
	    'ytick',[],...
	    'Fontsize',14,...
	    'Color',[0 0 0]);

    if s_fig>0
      fgnm=sprintf('%sGLBb08_HTC_Z%4i_Trf%3.3i_%i',...
		   pthfig,abs(round(Zref)),Tref,year);
      fprintf('Saving %s\n',fgnm);
      print('-dpng','-r250',fgnm);
    end
  end
end

% Area-mean over specified regions:
% Calculate FWC budget 
% for regions
HCR = sub_define_regions(HH);
NR=length(HCR);
for iyr=1:nyrs
  year=YR(iyr);
  htc=HTC(iyr).HeatCnt_J_m2*cff-htc0;
  for jrg=1:NR
    IN=HCR(jrg).IN;
    htc_tot=nansum(nansum(htc(IN).*Acell(IN))); %GJ
    atot=sum(Acell(IN));
    HCR(jrg).Heat_area(iyr,1)=htc_tot/atot; % GJ/m2
  end
end

figure(31); clf;
%POS(1,:)=[0.08 0.7 0.8 0.25];
%POS(2,:)=[0.08 0.37 0.8 0.25];
%POS(3,:)=[0.08 0.05 0.8 0.25];
POS(1,:)=[0.07 0.6 0.42 0.32];
POS(2,:)=[0.07 0.15 0.42 0.32];
POS(3,:)=[0.57 0.6 0.42 0.32];
POS(4,:)=[0.57 0.15 0.42 0.32];
for jrg=1:NR
  nm=sprintf('Anom. Ht.Cont*%5.1d J/m2, Z=%i, Trf=%4.1f, %s',...
	     cff,Zref,Tref,HCR(jrg).Name);
  ps=POS(jrg,:);
  axes('Position',ps);
  hcnt=HCR(jrg).Heat_area;
%  plot(YR,hcnt,'linewidth',2);
  hb = bar(YR,hcnt,0.9);
  set(hb,'Facecolor',[0.9 0.4 0]);
  title(nm);
  yl1=floor(min(hcnt)-0.1*abs(min(hcnt)));
  yl2=ceil(max(hcnt)+0.1*max(hcnt));
%  yl1=-0.63;
%  yl2=0.63;
  set(gca,'xlim',[1995.5 2015.5],...
	  'xtick',[1993:1:2015],...
	  'ylim',[yl1 yl2],...
	  'ytick',[-0.7:0.1:0.7],...
	  'tickdir','out',...
	  'fontsize',14);
  
  title(nm);
  set(gca,'xlim',[1995.5 2015.5],...
	  'xtick',[1993:2:2015],...
	  'xminortick','on',...
	  'tickdir','out',...
	  'Fontsize',14);
end
bottom_text(txtb);

if s_fig>0
  fgnm=sprintf('%sGLBb08_AreaMnHeatCont_Z%4.4i',...
	       pthfig,abs(round(Zref)));
  fprintf('Saving %s\n',fgnm);
  print('-dpng','-r250',fgnm);
  print('-depsc2',fgnm);
end

% Plot regions
f_reg=0;
if f_reg>0
% Plot regions
  nf=32;
  nrplt=[];
  sub_plot_regions(nf,HFLUX,HH,LON,LAT,nrplt);
  bottom_text(txtb);
  
  if s_fig>0
  fgnm=sprintf('%sGLBb08_AreaMnHeatCont_RegionsSubPol',...
	       pthfig);
    fprintf('Saving %s\n',fgnm);
    print('-dpng','-r250',fgnm);
  end

end










