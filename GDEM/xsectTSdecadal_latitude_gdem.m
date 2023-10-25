% Plot climatology T/S vertical sections
%
% Read instant output archive files:
% Note that original global HYCOM GLBb0.08 outputs 
% have been "remapped" onto ARCc0.08 grid
% see ../REMAP_ARCc/remap_gridGLB2ARC_archv
% Plotting was modified by L. Stefanova based on my codes
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

s_mat  = 0; % overridden if  s_extr=0
s_fig = 0;

pthbin  = '/Net/mars/ddmitry/hycom/GLBb2ARCc0.08/mean/'; % INPUT data GLBb0.08
pthmat  = '/nexsan/people/takis/lydia/HYCOM/data_mat/'; 
pthmat2 = '/Net/ocean/ddmitry/HYCOM/ARCc/data_mat/';
pthfig  = '/Net/ocean/ddmitry/HYCOM/ARCc/fig_GDEM4/'
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/'; 
pthin   = '/Net/data/GDEM4/';  % climatology data with no land

fmat    = sprintf('%sADPTH_decadal.mat',pthmat); 

ftopo = sprintf('%sdepth_ARCc0.08_07.nc',pthtopo); 
HH  = nc_varget(ftopo,'Bathymetry'); 
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(HH); 

% Convert longitudes to -180 +180 rrange
% piece inserted from xsectAO.m

elon=LON;alat=LAT;

A=elon;
[my,nx]=size(elon);
for i=1:nx
for j=1:my
  long=elon(j,i);
  dmm=long;
  while (abs(dmm)>180)
    if dmm<0
      dmm=dmm+360;
    else
      dmm=dmm-360;
    end;
  end;  % while dmm
  elon(j,i)=dmm;
end;  % for j
end;  % for i
clear A;
% Grid cell spacing
[DX,DY]=sub_dx_dy(LON,LAT);

xname = 'BeaufKara';

switch(xname);
  case('NAtl_Lat647');
   yS = 64.67;
   xS = -13;
   yE = 64.67;
   xE = 10;
   id='64.67N_13W-10E'

% Bering-to-Beaufort B
  case('BeaufBering');
   yS=72.6;
   xS=-128;
   yE=69.;
   xE=-172.5;
   id='72.6N-128W_69.0N_172.5W';

% Bering-to-Beaufort A
%yS=69.;
%xS=-172.5;
%yE=80.;
%xE=175;
%id='69N-172.5W_80N_175E';

 case('BeaufKara')
   yS=70;
   xS=-150;
   yE=73;
   xE=70;
   id='BeaufKara';
end  

dd = distance_spheric_coord(yS,xS,LAT,LON);
[jS,iS] = find(dd==min(min(dd)),1);
dd = distance_spheric_coord(yE,xE,LAT,LON);
[jE,iE] = find(dd==min(min(dd)),1);

% Boundary

%jS=25; iS=248; jE=25; iE=1016;
%id='_boundary'

fmat    = sprintf('%sXsect%s.mat',pthmat,id); % LS define the file for storing mat output


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
IJ=IJs;




%-------------------------------
   im=1
   finT=sprintf('%sptgdemv4f%2.2i.nc4',pthin,im);  % ptgdemv4f##.nc4
   zz=nc_varget(finT,'Depth');


  % loading data
  if s_mat==0
     fmatT=sprintf('%sTave_GDEM_00_full.mat',pthmat)
     fmatS=sprintf('%sSave_GDEM_00_full.mat',pthmat)
     %fmatDP=sprintf('%sDPave_GDEM_00_full.mat',pthmat)
     load(fmatT,'Tav');
     load(fmatS,'Sav');
     %load(fmatDP,'DPav');
  end;

% calculate depths from thicknessess

  %[ZZav,ZMav] = sub_thck2dpth(DPav);

  

% Depth along transsect
  for cntr=1:length(IJ)
    ii=IJ(cntr,1);
    jj=IJ(cntr,2);
     deepest(cntr)=min(find(-zz<HH(jj,ii))); 
  end;

  for cntr=1:length(IJ)
    ii=IJ(cntr,1);
    jj=IJ(cntr,2);

    Temperature(:,cntr)=squeeze(Tav(:,jj,ii)); 
    Salinity(:,cntr)=squeeze(Sav(:,jj,ii)); 
    DepthZ(:,cntr)=-zz; 
    
    Salinity(deepest(cntr):end,cntr)=nan;
    Temperature(deepest(cntr):end,cntr)=nan;
    DepthZ(deepest(cntr):end,cntr)=nan;

  end;
     %    Density=sw_dens0(Salinity,Temperature);
     %    Fld='Rho';c1=1025.1;c2=1028.1;
     %    sub_plot_xsectTSdecadal_gdem(Density,DepthZ,elon,alat,Fld,IJ,c1,c2,pthfig,id)
% plot temperature
Fld='T';
c1=-0.6;
c2=1.4;
% Fld='T';c1=0;c2=15;
sub_plot_xsectTSdecadal_gdem(Temperature,DepthZ,...
	elon,alat,Fld,IJ,c1,c2,pthfig,id,s_fig)

% plot salinity
Fld='S'; 
c1=34.4; 
c2=35;
% Fld='S';c1=34.5;c2=36;
sub_plot_xsectTSdecadal_gdem(Salinity,DepthZ,...
       elon,alat,Fld,IJ,c1,c2,pthfig,id,s_fig)

  


