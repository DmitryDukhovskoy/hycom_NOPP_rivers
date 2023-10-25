% Arctic Ocean 0.04 HYCOM-CICEv5 GOFS3.5
% Greenland runoff, no passive tracers
% analyze monthly mean sea ice fields
% extract from daily instantenous output
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

expt = 22;
pthout = '/Net/kronos/ddmitry/hycom/ARCc0.04/datamat/cice_mnth/';
pthmat = pthout;

f_conc = 0;
f_uv   = 0;

mmean = 1; % = 0 - instanteneous or monthly mean fields
yr   = 2017;
mo   = 8;
mday = 15;
YR   = yr;


pthbin = sprintf('/nexsan/people/ddmitry/hycom/ARCc0.04_%3.3i/data/%i_cice/',expt,YR);
%pthbin = '/nexsan/people/ddmitry/hycom/ARCc0.04_022/data_test/';  % test sim with corrected ice-ocean coupling
%pthbin = '/Net/kronos/ddmitry/hycom/ARCc0.04/cice_0901gofs35/'; % sea ice restart from GOFS3.5 GLBc
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';

btx = 'plot_daily_seaice.m';


ftopo = sprintf('%sdepth_ARCc0.04_17DD.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
Lmsk = HH*0;
Lmsk(HH<0)=1;


dnmb=datenum(YR,mo,mday);
iday=dnmb-datenum(YR,1,1)+1;

if mmean 
  flnm = sprintf('%s022_cice.%i-%2.2i-%2.2i.nc',pthbin,YR,mo,mday);
else
  flnm = sprintf('%s022_cice_inst.%i-%2.2i-%2.2i-00000.nc',pthbin,YR,mo,mday);
end

fprintf('Reading %s \n',flnm);
Hi = squeeze(nc_varget(flnm,'hi'));
Si = squeeze(nc_varget(flnm,'hs'));
Ci = squeeze(nc_varget(flnm,'aice'));
%Ui = squeeze(nc_varget(flnm,'uvel'));
%Vi = squeeze(nc_varget(flnm,'vvel'));
Hi(Hi==0)=nan;
Ci(HH>=0)=nan;
TLON = nc_varget(flnm,'TLON');
TLAT = nc_varget(flnm,'TLAT');


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
fprintf('Plotting H ice\n');
sub_plot_hice(Hi,Lmsk,xl1,xl2,yl1,yl2,c1,c2,fn,cmp);
contour(Ci,[0.15 0.15],'r','Color',[0.8 0.3 0],'Linewidth',1.5);

%stl = sprintf('CPL-TEST HYCOM-CICEv5, Hice, Cice=15%%, %2.2i/%2.2i/%4.4i',mday,mo,YR);
%stl = sprintf('HYCOM-CICEv5, Hice, Cice=15%%, %2.2i/%2.2i/%4.4i',mday,mo,YR);
stl = sprintf('GLBc CICE rest 0901, HYCOM-CICEv5, Hice, Cice=15%%, %2.2i/%2.2i/%4.4i',mday,mo,YR);
title(stl);
bottom_text(btx,'pwd',1);

if f_conc==1
%
% Sea ice concentration
%
cc1=0;
cc2=1;
pos=[0.12 0.1 0.82 0.82];
%CMP = colormap_PBYR(200,c1,c2);
CMP = colormap_blue_cyan_white(200,cc1,cc2);
cmp = CMP.colormap;
cmp = flipud(cmp);


fn=3;
fprintf('Plotting C ice\n');
dCi=Ci;
dCi(Ci<1e-9)=nan;
sub_plot_hice(Ci,Lmsk,xl1,xl2,yl1,yl2,cc1,cc2,fn,cmp);
contour(Ci,[0.15 0.15],'r','Color',[0 0 0],'Linewidth',1.5);
contour(Ci,[0.5 0.5],'r','Color',[0 0 0],'Linewidth',1.);

%stl = sprintf('CPL-TEST HYCOM-CICEv5, Hice, Cice=15%%, %2.2i/%2.2i/%4.4i',mday,mo,YR);
stl = sprintf('HYCOM-CICEv5, Cice, %2.2i/%2.2i/%4.4i',mday,mo,YR);
title(stl);
bottom_text(btx,'pwd',1);
end


if f_uv==1
%
% U ice
%
cs1=0;
cs2=0.4;
CMP = create_colormap_WBYR(200,c1,c2);
cmu = CMP.colormap;
Spi = sqrt(Ui.^2+Vi.^2);


fn=2;
fprintf('Plotting Ice U,V\n');
sub_plot_uice(Spi,Lmsk,xl1,xl2,yl1,yl2,cs1,cs2,fn,cmu);
%contour(Ci,[0.15 0.15],'r','Color',[0.8 0.3 0],'Linewidth',1.5);
% Plot ice u,v vectors
II=find(Ci<1e-3);
Ui(II)=nan;
Vi(II)=nan;
[mm,nn]=size(Ui);


scl=200;
cf=0.3;
beta=20;
lwd=1.0;
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
stl = sprintf('HYCOM-CICEv5, Uice, %2.2i/%2.2i/%4.4i',mday,mo,YR);
title(stl);

bottom_text(btx,'pwd',1,'fontsize',8);

end









