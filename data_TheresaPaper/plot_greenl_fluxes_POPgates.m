% Plot fluxes on the Greenland shelf
% calculated in greenl_fluxes_POPgates
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/Wavelet;
startup;

close all
clear


expt=110;
TV=11;
YR1=2008;
YR2=2008;


Cp = 4200; % J/kg K
Tref1= -1.8; % Ref T to calc. H flux
Tref2= 0; % Ref T to calc. H flux
hgg=1e20;
f_zgrd=0;  % =1 - calculate fluxes from z-grid interpolated U,T,S - less accurate
           % mostly for comparison and validation

pthmat = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_theresa/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';

btx='greenl_fluxes_POPgates.m';

fprintf('arc08-%3.3i Heat and Vol fluxes Greenland Shelf gates %i-%i\n',...
	expt,YR1,YR2);


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

% Vol flux for box on a eastern Gr Shelf
% Tingmiarmiut and Denmark Str. and 800m isobath
vFlx1=[];
vFlx2=[];
GvFlx=[];
for YR=YR1:YR2
  fmatout=sprintf('%shycom008_%3.3i_Greenl_flx_POPgates_%4.4i.mat',...
		    pthmat,expt,YR);
  fprintf('Loading %s\n',fmatout);
  load(fmatout);
  
  vFlx1=[vFlx1,VHFLX(3).VolFlxGrSh_m3s];
  vFlx2=[vFlx2,VHFLX(4).VolFlxGrSh_m3s];
%
% Flux across the Greenland contour:
  fmat2=sprintf('%s%3.3i_GreenlCntr_HVflx_daily_%4.4i.mat',...
		pthmat,expt,YR);
  fprintf('Loading %s\n',fmat2);
  load(fmat2);
  
  dmm=HFLX.VolFlx;
  igr1=VHFLX(3).GrCntr_I;
  jgr1=VHFLX(3).GrCntr_J;
  igr2=VHFLX(4).GrCntr_I;
  jgr2=VHFLX(4).GrCntr_J;
  
  if ~exist('Igr1','var') 
    D=sqrt((Ig-igr1).^2+(Jg-jgr1).^2);
    Igr1=find(D==min(D));
    D=sqrt((Ig-igr2).^2+(Jg-jgr2).^2);
    Igr2=find(D==min(D));
  end  
    
  vv1=nansum(dmm(:,Igr1:Igr2),2);
  GvFlx=[GvFlx,vv1'];
  
end

Gr=GvFlx;
Gr(abs(Gr)>1.8e7)=mean(Gr);
nn=length(Gr);
[B,A]=butter(5,60/nn,'low');
Grf=filtfilt(B,A,Gr);
vf1=filtfilt(B,A,vFlx1);
vf2=filtfilt(B,A,vFlx2);

figure(1); clf;
plot(vf1); 
hold on;
plot(-vf2); % the flux is + into the box, N gate
plot(Grf);
plot(vf2-vf1); % what is needed to compensate the difference btw gate fluxes

legend('SGate','NGate','Contour','-(S+N)');
title('Filtered Vol Fluxes, Denmark box');
