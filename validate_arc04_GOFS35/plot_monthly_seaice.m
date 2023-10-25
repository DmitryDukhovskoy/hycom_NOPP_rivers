% Plot monhtly CICE output 
% extracted in monthly_seaice.m
%
% Arctic Ocean 0.04 HYCOM-CICEv5 GOFS3.5
% Greenland runoff, no passive tracers
% analyze monthly mean sea ice fields
% extract from daily instantenous output
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

expt = 22;
pthout = '/Net/kronos/ddmitry/hycom/ARCc0.04/datamat/cice_mnth/';
pthmat = pthout;
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';

btx = 'plot_monthly_seaice.m';


ftopo = sprintf('%sdepth_ARCc0.04_17DD.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
Lmsk = HH*0;
Lmsk(HH<0)=1;

[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2

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


% Plot maps:
yr1=2016;
yr2=2017;
icc=0;
for YR=yr1:yr2
  pthbin = sprintf('/nexsan/people/ddmitry/hycom/ARCc0.04_%3.3i/data/%i_cice/',expt,YR);

  for im=1:3:12
    fmatout = sprintf('%scice_monthly_%i%2.2i.mat',pthmat,YR,im);
    if ~exist(fmatout), continue; end;
    fprintf('Loading %s\n',fmatout);
    load(fmatout);

    icc=icc+1;

    ci = ICEMO.Ice_conc;
    si = ICEMO.Snow_thck;
    hi = ICEMO.Ice_thck;
    hi(hi==0)=nan;
    ci(HH>=0)=nan;
    
%    hi(ci<0.15)=nan;  % cutoff below 15% concentration
    fn=icc;
    sub_plot_hice(hi,Lmsk,xl1,xl2,yl1,yl2,c1,c2,fn,cmp);
    contour(ci,[0.15 0.15],'r','Color',[0.8 0.3 0],'Linewidth',1.5);
%    contour(si,[0.05 0.05],'k-','Color',[0.95 0.95 0.95]);
%    contour(si,[0.1 0.1],'k-','Color',[1 1 1]);
%    contour(si,[0.15 0.15],'k-','Color',[1 1 1],'Linewidth',1.5);
%    contour(si,[0.2 0.2],'k-','Color',[1 1 1],'Linewidth',2);
%    contour(si,[0.25 0.25],'k-','Color',[1 1 1],'Linewidth',1.5);
  
    stl = sprintf('HYCOM-CICEv5, Hice, Cice=15%%, %2.2i/%4.4i',im,YR);
    title(stl);
    bottom_text(btx,'pwd',1);

%    I = find(~isnan(hi));
%    hsnow = si(I);
%    hice  = hi(I);

%    keyboard

  end
end




