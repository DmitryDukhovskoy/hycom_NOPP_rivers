% CICE5 is mapd from GLBc0.04 GOFS3.5 CICE5
% using D. Hebert code
% on Gordon
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear


pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/topo_grid/';
pthdat = '/Net/kronos/ddmitry/hycom/ARCc0.04/cice/';

dnmb = datenum(2017,09,01,9,0,0);
dv = datevec(dnmb);
%frst = sprintf('%scice.restart.2017010109_interp.nc',pthdat);
frst = sprintf('%scice.restart.%4.4i%2.2i%2.2i%2.2i_interp.nc',pthdat,dv(1:4));

ftopo = sprintf('%s/depth_ARCc0.04_17DD.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
Lmsk = HH*0;
Lmsk(HH<0)=1;

%
%uu = nc_varget(frst,'uvel');

% Plot conc:
ncat = 5;
dmm = nc_varget(frst,'aicen');
aice = squeeze(sum(dmm,1));
clear dmm

btx = 'plot_restart_cice5.m';


% Ice thickness:
% volume per 1 m2:
dmm = nc_varget(frst,'vicen');
hice = squeeze(sum(dmm,1));
hice(hice==0)=nan;
clear dmm

%cmp = colormap('jet');

figure(1);
nf=1;
c1=0;
c2=4;
pos=[0.12 0.1 0.82 0.82];
xl1= 50;
xl2= nn;
yl1= 500;
yl2= 4000;
%CMP = colormap_PBYR(200,c1,c2);
CMP = create_colormapBGY(200,c1,c2);
cmp = CMP.colormap;


fn=1;
sub_plot_hice(hice,Lmsk,xl1,xl2,yl1,yl2,c1,c2,fn,cmp);
contour(aice,[0.15 0.15],'r','Color',[0.8 0.3 0],'Linewidth',1.5);

mday = dv(3);
mo = dv(2);
YR = dv(1);
%stl = sprintf('CPL-TEST HYCOM-CICEv5, Hice, Cice=15%%, %2.2i/%2.2i/%4.4i',mday,mo,YR);
stl = sprintf('GLBc GOFS3.5 restart for ARCc0.04 Hice, Cice=15%%, %2.2i/%2.2i/%4.4i',mday,mo,YR);
title(stl);
bottom_text(btx,'pwd',1);


%
%    Ice Concentration
%
cc1=0;
cc2=1;
pos=[0.12 0.1 0.82 0.82];
%CMP = colormap_PBYR(200,c1,c2);
CMP = colormap_blue_cyan_white(200,cc1,cc2);
cmp = CMP.colormap;
cmp = flipud(cmp);


fn=3;
Ci = aice;
Ci(HH>=0) = nan;
sub_plot_hice(Ci,Lmsk,xl1,xl2,yl1,yl2,cc1,cc2,fn,cmp);
contour(Ci,[0.15 0.15],'r','Color',[0 0 0],'Linewidth',1.5);
contour(Ci,[0.5 0.5],'r','Color',[0 0 0],'Linewidth',1.);

%stl = sprintf('CPL-TEST HYCOM-CICEv5, Hice, Cice=15%%, %2.2i/%2.2i/%4.4i',mday,mo,YR);
stl = sprintf('GLBc GOFS3.5 restart for 0.04ARCc, Cice, %2.2i/%2.2i/%4.4i',mday,mo,YR);
title(stl);
bottom_text(btx,'pwd',1);



%
% U ice
%
Ui = nc_varget(frst,'uvel');
Vi = nc_varget(frst,'vvel');
cs1=0;
cs2=0.8;
CMP = create_colormap_WBYR(200,c1,c2);
cmu = CMP.colormap;

II=find(aice<1e-3);
Ui(II)=nan;
Vi(II)=nan;
JJ=find(HH>=0);
Ui(JJ)=nan;
Vi(JJ)=nan;
Spi = sqrt(Ui.^2+Vi.^2);

[mm,nn]=size(Ui);

fn=2;
sub_plot_uice(Spi,Lmsk,xl1,xl2,yl1,yl2,cs1,cs2,fn,cmu);
%contour(Ci,[0.15 0.15],'r','Color',[0.8 0.3 0],'Linewidth',1.5);
% Plot ice u,v vectors
scl=100;
cf=0.3;
beta=20;
lwd=1.2;
v_col=[0 0 0];
dii=40;
for ii=dii:dii:nn
  for jj=dii:dii:mm
    clear u v
    u = Ui(jj,ii);
    v = Vi(jj,ii);
    if isnan(u), continue; end;

    sp = sqrt(u*u+v*v);


    x0=ii;
    y0=jj;

    x1=x0+u*scl;
    y1=y0+v*scl;
    draw_arrowF(x0,x1,y0,y1,cf,beta,v_col,lwd);

  end
end
%stl = sprintf('CPL-TEST HYCOM-CICEv5, Uice, %2.2i/%2.2i/%4.4i',mday,mo,YR);
stl = sprintf('GLBc GOFS3.5 restart for ARCc0.04 Uice %2.2i/%2.2i/%4.4i',mday,mo,YR);
title(stl);
bottom_text(btx,'pwd',1);










