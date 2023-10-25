% Plot FW content - calculated in fwc_natl.m + fwc_natl2016.m
% relative to Sref
% Using same regions as in EN4 analysis
% integrating to the depth of Sref
% Following Haine
% across straits 
% /nexsan/GLBa0.08/expt_90.9/data
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

s_fig  = 0;

rg=9806;  % convert pressure to depth, m
Sref=35; % N.Atl. is too saline
TV = '07';

pth190='/nexsan/GLBb0.08/GLBb0.08_190/data/meanstd/'; % 1993-1994
pth191='/nexsan/GLBb0.08/GLBb0.08_191/data/meanstd/'; % 1995-2012
pthmat = '/Net/mars/ddmitry/hycom/GLBb0.08/data_mat/';
pthfig = '/Net/mars/ddmitry/hycom/GLBb0.08/fig_fw/';
ftopo = sprintf('%sGLBb_T07_subset_Natl.mat',pthmat);
pthglb = '/nexsan/GLBb0.08/GLBb0.08_191/topo/';
%fmat = sprintf('%sFWC_natl_1993-2015.mat',pthmat);
fmat = sprintf('%sFWC_natl_1993-2016.mat',pthmat);
txtb = 'plot_fwc_natl_v2.m';

fprintf('Loading topo %s\n',ftopo);
load(ftopo);

% Grid cell spacing
[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2

[mm,nn]=size(HH);
%[X,Y]=meshgrid((1:nn),(1:mm));


fprintf('Loading %s\n',fmat);
load(fmat);
nyrs=length(FWC);
for iy=1:nyrs
  YR(iy,1)=FWC(iy).Year;
end


% Calculate time-mean:
fwc0=zeros(mm,nn);
cc=0;
for iy=1:7
  cc=cc+1;
  fwc=FWC(iy).Fwc_m; % delta htc
  fwc(isnan(fwc))=0;
  fwc(HH>=0)=nan;
  fwc0=fwc0+fwc;
end
fwc0=fwc0/cc;

nint=200;
c1=0;
c2=20;
CMP = create_colormap_freshwater(nint,c1,c2);
cnt=CMP.intervals;
cmp=CMP.colormap;

%YRPLT=[1993:2:2015];
YRPLT=[];
npl=length(YRPLT);
cc=0;
if npl>0
  for iyr=1:npl
    year=YRPLT(iyr);
    iy=find(YR==year);
    fwc=FWC(iy).Fwc_m;
    fwc(fwc<1e-1)=nan;
%    fwc=fwc-fwc0; % anomaly

%    cc=cc+1;
    figure(1); clf;
    hold on;
    contour(HH,[0 0],'k','Linewidth',1.6);
    contour(HH,[-5000:1000:-20],'Color',[0.8 0.8 0.8],'Linewidth',1.);
    pcolor(fwc); shading flat;
    colormap(cmp);
    caxis([c1 c2])
%    colorbar
    ttl=sprintf('GLBb0.08 FWC, m, Sref %4.1f, %i',Sref,year);
    title(ttl);
    axis('equal');
    set(gca,'xlim',[1 nn],...
	'ylim',[1 mm],...
	'xtick',[],...
	'ytick',[]);
    txtb='plot_fwc_natl_v2.m';
    bottom_text(txtb,'pwd',1);
    
  hght=[];
  lngth=[];
  mint=20;
  mbx=mint;
  fsz=13;
  bxc='k';
  posc=[0.78 0.11 0.8 0.08];
  aend=0;
  [az,axc]  = colorbar_vert (cmp,cnt,hght,lngth,mint,fsz,bxc,posc,mbx,aend);

    if s_fig>0
      fgnm=sprintf('%sGLBb08_FWC_Srf%3.3i_%i',...
		   pthfig,Sref,year);
      fprintf('Saving %s\n',fgnm);
      print('-dpng','-r250',fgnm);
    end
  end
end



% Calculate FWC budget 
% for regions
%keyboard
FWB = sub_define_NatlRegions(HH,LON,LAT);
NR = length(FWB);
for iyr=1:nyrs
  year=YR(iyr);
  fwc=FWC(iyr).Fwc_m;
  for jrg=1:NR
    IN=FWB(jrg).IN;
    fwvol=nansum(nansum(fwc(IN).*Acell(IN)*1e-9)); %km3
    area = sum(Acell(IN))*1e-6; % km2
    FWB(jrg).FWvol_km3(iyr,1)=fwvol; % km3
    FWB(jrg).FWvol_m(iyr,1)=fwvol/area*1e3; % m of Fwater
  end
end

xlbl=[];
cll=0;
for yr=YR(1):YR(end);
  cll=cll+1;
  if mod(yr,5)==0
    xlbl{cll}=sprintf('%i',yr);
  else
    xlbl{cll}=' ';
  end
end

% Plot annual FWC
for jrg = 1:NR
  figure(jrg); clf;
  axes('position',[0.1 0.5 0.85 0.4]);
  nm=sprintf('FWC, m, Sref=%4.1f %s',Sref,FWB(jrg).Name);
  fwv=FWB(jrg).FWvol_m;
  fwv=FWB(jrg).FWvol_km3;
  hb=bar(YR,fwv,0.9);
  set(hb,'Facecolor',[0.3 0.7 1]);
  title(nm,'Fontsize',11);
  yl1=floor(min(fwv)-0.1*min(fwv));
  yl2=ceil(max(fwv)+0.1*max(fwv));
  dyl=1;
  if (yl2-yl1)>10,
    dyl=2;
  end
  
  set(gca,'xlim',[YR(1)-0.5 YR(end)+0.5],...
	  'xtick',[YR(1):YR(end)],...
	  'xticklabel',xlbl,...
	  'ylim',[yl1 yl2],...
	  'ytick',[0:dyl:50],...
	  'tickdir','out',...
	  'fontsize',14);
%  ylabel('FWC, km3');

  bottom_text(txtb,'pwd',1,'position',[0.08 0.3 0.7 0.2]);

  if s_fig>0
    fgnm=sprintf('%sGLBb08_annFWC_%s_Srf%3.3i',...
		 pthfig,FWB(jrg).Name,Sref);
    fprintf('Saving %s\n',fgnm);
    print('-dpng','-r250',fgnm);
  %  print('-depsc2',fgnm);
  end

end

% Plot regions
f_reg=0;
if f_reg>0
% Plot regions
  nf=10;
  nrplt=[];
  sub_plot_regions(nf,FWB,HH,LON,LAT,nrplt);
  bottom_text(txtb);
  
  if s_fig>0
    fgnm=sprintf('%sAOFlx_regions',pthfig);
    fprintf('Saving %s\n',fgnm);
    print('-dpng','-r250',fgnm);
  end

end




