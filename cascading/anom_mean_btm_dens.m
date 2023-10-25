% Anomalies and mean bottom densities
% monthly mean near-bottom densities
% are calculated in bottom_dens.m
% 
% expt_110 - no Greenland runoff
% expt 112 - with Greenland runoff
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

yr1=2005;
yr2=2009;

regn = 'ARCc0.08';
%expt = 110; % no Greenland runoff  
expt = 112;  % Greenland runoff

rg = 9806;
% Thershold dRho
dRho0=1e-4*1023;


pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat2/',expt);

ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[XX,YY]=meshgrid((1:nn),(1:mm));


fprintf('expt %3.3i\n\n',expt);

frho=sprintf('%srho_btm_mean.mat',pthmat);
f_mean=0;
if f_mean==1
  for yr=yr1:yr2
    fmat = sprintf('%sarc08_%3.3i_btmRho_%i.mat',pthmat,expt,yr);
    fprintf('Load %s\n',fmat);
    load(fmat);

    cntr=HH*0;
    rho_sum=HH*0;
    for mo=1:12
      rho=RBTM(mo).Rho;
      I=find(rho>999.9);
      rho_sum(I)=rho_sum(I)+rho(I);
      cntr(I)=cntr(I)+1;
    end
  end

  Rho_mean=rho_sum./cntr;

  frho=sprintf('%srho_btm_mean.mat',pthmat);
  fprintf('Saving %s\n',frho);
  save(frho,'Rho_mean');
else
  fprintf('Loading mean rho: %s\n',frho);
  load(frho);
end

MM={'J','F','M','A','M','J','J','A','S','O','N','D'};

sM=[];

%keyboard;
cc=0;
for yr=yr1:yr2
  fmat = sprintf('%sarc08_%3.3i_btmRho_%i.mat',pthmat,expt,yr);
  fprintf('Load %s\n',fmat);
  load(fmat);

  cntr=HH*0;
  rho_sum=HH*0;
  for mo=3:3
    rho=RBTM(mo).Rho;
    I=find(rho<999.9);
    rho(I)=nan;
    anom_rho=rho-Rho_mean;
    I=find(~isnan(rho));
    rho_sum(I)=rho_sum(I)+anom_rho(I);
    cntr(I)=cntr(I)+1;

    cc=cc+1;
    if yr==yr1
     sM=[sM,MM{mo}];
    end
  end
end

arho=rho_sum./cntr;

% Colormap:
c1=0;
c2=0.3;
nint= 200;
CMP = create_colormap5(nint,c1,c2);
cmp0= CMP.colormap;
for k=1:20
  cmp0(k,:)=[1 1 1];
end
nav = 15;
cmp = smooth_colormap(cmp0,nav);
%cmp=cmp0;
cnt = CMP.intervals;

hmsk=HH;
hmsk(HH<0)=nan;

figure(1); clf;
pcolor(arho); shading flat;
colormap(cmp);
caxis([c1 c2]);

axis('equal');
set(gca,'xlim',[450 1000],...
        'ylim',[380 1100],...
	'xtick',[],...
	'ytick',[]);

hold on;
contour(HH,[-800 -800],'Color',[0.4 0.4 0.4]);

freezeColors;
pcolor(hmsk); shading flat;
colormap([0.8 0.8 0.8]);
freezeColors;

colormap(cmp);
hb=colorbar;
set(hb,'Position',[0.9 0.1 0.02 0.8],...
       'Ticklength',0.02,...
       'Fontsize',14);

stl=sprintf('HYCOM0.08-%3.3i, dRho_bottom, %i-%i, %s',...
             expt,yr1,yr2,sM);
title(stl,'Interpreter','none'); 

btx='anom_mean_btm_dens.m';
bottom_text(btx,'pwd',1);









  
  
