% Average heat and volume fluxes across Greenland contour
% extracted in ocn_hflx_greenl008.m

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/seawater;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

YR1=1993;
YR2=2009;

regn = 'ARCc0.08';
expt = 110; % experiment without runoff
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
%pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
pthmat  = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_mat2/';
btx='HtVlm_flux_contour_mean_season.m';

ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);

GC = sub_greenl_isobath(HH,LON,LAT);
Hs = GC.Hbottom; % bottom along section

% ================================================
% Plot Greenland map and the contour
% ================================================
f_pltgr=0;
if f_pltgr==1
  fn=10;
  sub_plot_Greenl_contour(HH,LON,LAT,fn,GC);

  bottom_text(btx,'pwd',1);
end

% Filter spatially along the contour
% to get rid of highly oscillatory flux
Wn = 1/40;
[Bf,Af] = butter(9,Wn,'low');

ntot=0;
ndjf=0;
njja=0;
nmam=0;
nson=0;
for iyr=YR1:YR2
  fmat=sprintf('%s%3.3i_Greenl_HVflx_%i.mat',...
                   pthmat,expt,iyr);
  fprintf('Loading %s\n',fmat);
  load(fmat);
  if ntot==0; 
    Vtot=HFLX(1).Vol_flux_m3s*0; 
    Htot=HFLX(1).Vol_flux_m3s*0; 
    Vdjf=Vtot;
    Vmam=Vtot;
    Vjja=Vtot;
    Vson=Vtot;
    Hdjf=Vtot;
    Hmam=Vtot;
    Hjja=Vtot;
    Hson=Vtot;
    nsgm=length(Vtot);
    DL=HFLX(1).DistCntr; % distances, m
    dX=diff(DL); % segment lengths, m
    dX(nsgm)=dX(nsgm-1); 
    Hb=HFLX(1).Hbottom;
  end;
  
  for imo=1:12
    ntot=ntot+1;
% Correct overall transport to make=0    
    V=HFLX(imo).Vol_flux_m3s;
    dV=nansum(V);
    Asgm=abs(Hb).*DL; % m2
    Asgm(Hb>0)=nan;
    At=nansum(Asgm);
    dltV=-dV/At*Asgm;
    vtrt=(HFLX(imo).Vol_flux_m3s+dltV)./DL; % m3/s per 1 m
    
    Vtot=Vtot+vtrt;
    ZZ=HFLX(imo).ZZ;
    [a,b]=size(ZZ);
    ZZ=[zeros(1,b);ZZ];
    dZ=[];
    for k=1:a
      dZ(k,:)=abs(ZZ(k+1,:)-ZZ(k,:));
    end
%
% Depth integrated:
% Correct heat flux
% based on transport flux correction
    frV=dltV./abs(V); % fraction corrected
    frV(V==0)=nan;
    HF=HFLX(imo).Hflux_W;
    HFs=HF;
    for k=1:a-1
      HFs(k,:)=HF(k,:)+HF(k,:).*frV;
    end
    Hf=nansum((HFs./DL).*dZ,1); % W/m
    Htot=Htot+Hf;
%    
% Averaging by seasons
    if imo==12 | imo==1 | imo==2
      ndjf=ndjf+1;
      Vdjf=Vdjf+vtrt;
      Hdjf=Hdjf+Hf;
    end
    if imo>=3 & imo<=5
      nmam=nmam+1;
      Vmam=Vmam+vtrt;
      Hmam=Hmam+Hf;
    end
    if imo>=6 & imo<=8
      njja=njja+1;
      Vjja=Vjja+vtrt;
      Hjja=Hjja+Hf;
    end
    if imo>=9 & imo<=11
      nson=nson+1;
      Vson=Vson+vtrt;
      Hson=Hson+Hf;
    end
  end

  fprintf('%i done, ntot=%i, ndjf=%i,dmam=%i, njja=%i, nson=%i\n',...
	iyr, ntot, ndjf, nmam, njja, nson);
  
end

% Average:
Vtot=Vtot/ntot;
Htot=Htot/ntot;
Vdjf=Vdjf/ndjf;
Hdjf=Hdjf/ndjf;
Vmam=Vmam/nmam;
Hmam=Hmam/nmam;
Vjja=Vjja/njja;
Hjja=Hjja/njja;
Vson=Vson/nson;
Hson=Hson/nson;

% Filter along the contour:
f_flt=0;
if f_flt==1
  dmm = Vtot;
  dmm(isnan(dmm))=0;
  yy = filtfilt(Bf,Af,dmm); %
  yy(Hb>0)=nan;
  Vtot=yy;

  dmm = Htot;
  dmm(isnan(dmm))=0;
  yy = filtfilt(Bf,Af,dmm); %
  yy(Hb>0)=nan;
  Htot=yy;

  dmm = Vdjf;
  dmm(isnan(dmm))=0;
  yy = filtfilt(Bf,Af,dmm); %
  yy(Hb>0)=nan;
  Vdjf=yy;

  dmm = Hdjf;
  dmm(isnan(dmm))=0;
  yy = filtfilt(Bf,Af,dmm); %
  yy(Hb>0)=nan;
  Hdjf=yy;

  dmm = Vmam;
  dmm(isnan(dmm))=0;
  yy = filtfilt(Bf,Af,dmm); %
  yy(Hb>0)=nan;
  Vmam=yy;

  dmm = Hmam;
  dmm(isnan(dmm))=0;
  yy = filtfilt(Bf,Af,dmm); %
  yy(Hb>0)=nan;
  Hmam=yy;


  dmm = Vjja;
  dmm(isnan(dmm))=0;
  yy = filtfilt(Bf,Af,dmm); %
  yy(Hb>0)=nan;
  Vjja=yy;

  dmm = Hjja;
  dmm(isnan(dmm))=0;
  yy = filtfilt(Bf,Af,dmm); %
  yy(Hb>0)=nan;
  Hjja=yy;


  dmm = Vson;
  dmm(isnan(dmm))=0;
  yy = filtfilt(Bf,Af,dmm); %
  yy(Hb>0)=nan;
  Vson=yy;

  dmm = Hson;
  dmm(isnan(dmm))=0;
  yy = filtfilt(Bf,Af,dmm); %
  yy(Hb>0)=nan;
  Hson=yy;
end

% Get coordinates:
II=HFLX(1).GrCntr_II;
JJ=HFLX(1).GrCntr_JJ;
for ii=1:length(II);
  i=II(ii);
  j=JJ(ii);
  lon(ii)=LON(j,i);
  lat(ii)=LAT(j,i);
end

stl='Non filtered data';
if f_flt==1
  stl='Spatial filtering, Butterworth';
end

VHFLX.Title=sprintf('0.08 HYCOM-CICE-%3.3i %i-%i',expt,YR1,YR2);
VHFLX.Info=stl;
VHFLX.Hbottom=Hb;
VHFLX.Cntr_lon=lon;
VHFLX.Cntr_lat=lat;
VHFLX.Indx=II;
VHFLX.Jndx=JJ;
VHFLX.Segm_dx_m=DL; % segment length, m
VHFLX.Units='Transport: m3/s per 1 m, Heat flx: W/m';
VHFLX.Trt_tot=Vtot;
VHFLX.Hflx_tot=Htot;
VHFLX.Trt_djf=Vdjf;
VHFLX.Hflx_djf=Hdjf;
VHFLX.Trt_mam=Vmam;
VHFLX.Hflx_mam=Hmam;
VHFLX.Trt_jja=Vjja;
VHFLX.Hflx_jja=Hjja;
VHFLX.Trt_son=Vson;
VHFLX.Hflx_son=Hson;

fout=sprintf('%sarc08_%3.3i_Greenl_VolHeat_seasons_%i-%i.mat',pthmat,expt,YR1,YR2);
fprintf('Saving %s\n',fout);
save(fout,'VHFLX');

