addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

s_mat=0; % =0 - extract data, no save; =1 extract & save, =2 - load saved
s_fig=1;
rg=9806;  % convert pressure to depth, m
fld0='tracer'; %'salin' - background field
plr=20; % highlight this interface

pthbin  = '/Net/mars/ddmitry/hycom/GLBb2ARCc0.08/mean/';
pthmat  = '/Net/ocean/ddmitry/HYCOM/ARCc/GLBb0.08/data_mat/';
pthfig  = '/Net/ocean/ddmitry/HYCOM/ARCc/GLBb0.08/fig_climatology/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';

fmat    = sprintf('%sfwc_hycom.mat',pthmat);

ftopo = sprintf('%sdepth_ARCc0.08_07.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(HH);


% Grid cell spacing
[DX,DY]=sub_dx_dy(LON,LAT);

yS = 65.06;
xS = -20;
yE = 65.06;
xE = 20;

dd = distance_spheric_coord(yS,xS,LAT,LON);
[jS,iS] = find(dd==min(min(dd)),1);
dd = distance_spheric_coord(yE,xE,LAT,LON);
[jE,iE] = find(dd==min(min(dd)),1);

% Specify segments of the x-section
IJs=[iS, jS; ...
     iE, jE];
     
nij=size(IJs,1);
IIs=[];
JJs=[];
for ii=1:nij-1
  i1=IJs(ii,1);
  i2=IJs(ii+1,1);
  j1=IJs(ii,2);
  j2=IJs(ii+1,2);
  [I,J]=sub_xsct_indx(i1,j1,i2,j2);
  if size(I,1)==1;
    I=I';
    J=J';
  end
  
  IIs=[IIs;I];
  JJs=[JJs;J];
end;

IJs=[IIs,JJs];

nS=length(IIs);
clear Xl Yl
for ii=1:nS
  i0=IJs(ii,1);
  j0=IJs(ii,2);
  Xl(ii,1)=LON(j0,i0);
  Yl(ii,1)=LAT(j0,i0);
end;
INDs=sub2ind(size(HH),JJs,IIs);

f_map=0;
if f_map>0
  figure(1); clf;
  hold on
  contour(HH,[0 0],'k');
  contour(HH,[-5000:500:-100],'Color',[0.7 0.7 0.7]);
  plot(IIs,JJs,'b.-');
end






