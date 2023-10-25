% analyze fluxes extracted in fluxes_NorthAtlantic.m
% through N. Atlantic openings
% for budgets
%
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/Wavelet;
startup;

close all
clear


expt=112;
TV=11;
YR1=2007;
YR2=2007;

Sr1=34.8;
Sr2=34.95;
hgg=1e20;
f_zgrd=0;  % =1 - calculate fluxes from z-grid interpolated U,T,S - less accurate
           % mostly for comparison and validation

pthmat = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_natl/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';

btx='anls_fluxes_NorthAtlantic.m';


ftopo = sprintf('%sdepth_ARCc0.08_%2.2i.nc',pthtopo,TV); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);


FLX=[];
FLX(1).TM=[];
FLX(1).Vol=[];
FLX(1).FW1=[];
FLX(1).FW2=[];
for YR=2005:2007
  fmatout=sprintf('%shycom008_%3.3i_NAtl_fluxes_%4.4i.mat',...
		pthmat,expt,YR);
  fprintf('Loading %s\n',fmatout);
  load(fmatout,'VHFLX');
  
  for ik=1:7
    FLX(ik).Name=VHFLX(ik).Name;
    dmm=VHFLX(ik).VolFlx_m3s;
    dmm=dmm(:);
    bmm=FLX(ik).Vol;
    FLX(ik).Vol=[bmm;dmm];
    dmm=VHFLX(ik).Time;
    dmm=dmm(:);
    bmm=FLX(ik).TM;
    FLX(ik).TM=[bmm;dmm];
    dmm=VHFLX(ik).SFlx_S1;
    dmm=dmm(:);
    bmm=FLX(ik).FW1;
    FLX(ik).FW1=[bmm;dmm];
    dmm=VHFLX(ik).SFlx_S2;
    dmm=dmm(:);
    bmm=FLX(ik).FW2;
    FLX(ik).FW2=[bmm;dmm];
  end
    
end

fprintf(' --------------------------------- \n\n');

for ik=1:7
  nm=FLX(ik).Name;
  Vol=FLX(ik).Vol*1e-6;
  FW1=FLX(ik).FW1*1e-3;
  FW2=FLX(ik).FW2*1e-3;
  N=length(Vol);
  
  mvol=nanmean(Vol);
  p1=prctile(Vol,10);
  p9=prctile(Vol,90);
  verr=std(Vol)/sqrt(N);

  ms1=nanmean(FW1);
  p11=prctile(FW1,10);
  p19=prctile(FW1,90);
  s1err=std(FW1)/sqrt(N);
  
  ms2=nanmean(FW2);
  p21=prctile(FW2,10);
  p29=prctile(FW2,90);
  s2err=std(FW2)/sqrt(N);
  
  sv2km = 1e6*3600*24*365*1e-9; % Sv -> km3/yr
  msv2km = 1e3*3600*24*365*1e-9; % mSv -> km3/yr
  
  fprintf('%s: Vol=%5.2f+/-%5.2f Sv, %6.1f+/-%6.1f km3/y3\n',...
	  nm, mvol, verr, mvol*sv2km, verr*sv2km);
  fprintf('%s: S=%4.2f FWF=%5.2f+/-%5.2f mSv, %6.2f+/-%6.2f km3/y3\n',...
	  nm, Sr1,ms1, s1err, ms1*msv2km, s1err*msv2km);
  fprintf('%s: S=%4.2f FWF=%5.2f+/-%5.2f mSv, %6.2f+/-%6.2f km3/y3\n',...
	  nm, Sr2,ms2, s2err, ms2*msv2km, s2err*msv2km);
  fprintf(' --------------------------------- \n\n');
  
end
