% Updated code:
% Spatial map of dltS for specified day
% Uses corrected code extracting tracer
% that integrates tracer over the specified layers
% see: extr_MassTrcr_month
% Estimate S change in the upper layer due to surplus FW flux
% from Greenland 
% Note estimate of Greenland runoff fraction
% for specified boxes may use wrong
% code - tracer concentration averaged over layers
% should be integrated giving a total mass within
% a layer

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers
startup;

close all
clear

s_fig = 0;

dnmb = datenum(2016,1,30);
dv0  = datevec(dnmb);
yr   = dv0(1);
iday = dnmb-datenum(yr,1,1)+1;

nTr = 1; 
IntSrf = 1; % = 0 - integrates over all layers from 0 - ilv, for ilv=1 is the same
%ilv = 1; % 0-50m
%ilv = 2; % 50-150m
%ilv = 3; % 150-300m
ilv = 4; % 300-500 m
%ilv = 5; % whole depth <- do not use this

% Select figures to plot:
pfg1 = 0; % plot coefficient = Trac Mass in Layer/Total Trc Mass
pfg2 = 1; % plot FWC change rescaled to dzRef m for comparison
pfg3 = 1; % dlt S change within the layer


% Vertical layers
LRS = load('LRS.dat');
nlrs= length(LRS); 


zz1 = LRS(ilv,1);
zz2 = LRS(ilv,2);

if IntSrf==1; zz1=0; end;

dz=abs(zz2-zz1);

fprintf('S change, ilv=%i, %i - %i, nTr=%i, %i/%i/%i\n',ilv,zz1,zz2,nTr,dv0(1:3));

regn = 'ARCc0.08';
expt = 110;  
pthfig  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/fig_NAtl/',expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
pthsav  = sprintf('/Net/tholia/ddmitry/hycom/ARCc0.08/data_mat%3.3i/',expt);
%fmout   = sprintf('%sfraction_trcr_NAtl.mat',pthmat);

ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

% Grid cell spacing
[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2


hmsk=HH;
hmsk(HH<0)=nan;

%
% Get layer-averaged S - as a reference S
% to calculate dlt S
% prepared in anls_TS/mnthly_arc08_layers_S.m
YR=dv0(1);
mo=dv0(2);
fsout = sprintf('%sarc08_%3.3i_mnthS_lrs_%4.4i%2.2i.mat',pthsav,expt,YR,mo);
fprintf('Loading S averaged %s\n',fsout);
load(fsout);
Sav = meanS(ilv).Savrg;
Sav(Sav==0)=nan;
hZ = abs(LRS(ilv,2)-LRS(ilv,1));

if IntSrf==1
  smm = HH*0;
  shZ = 0; 
  for ikk=1:ilv
    Sav = meanS(ikk).Savrg;
    Sav(Sav==0)=nan;
    hZ = abs(LRS(ikk,2)-LRS(ikk,1));
    smm = smm+Sav.*hZ;
    shZ = shZ+hZ;
  end
  Sav = smm./shZ;
  hZ = dz;
end


%dSregions_map_v2


% Tracer fraction in grid cells
nTr = 1;
%rr = sub_fraction_tracer(expt,regn,Acell,HH,dnmb,nTr,ilv,hZ);
rr = sub_fraction_tracerMass(expt,regn,Acell,HH,dnmb,nTr,ilv,hZ,IntSrf);

frv = sprintf('%sGreenland_cumFWFlux.mat',pthmat);

f_griv = 0; %=1 - rederive cumulative Greenland FWFlux from Bamber
if f_griv==1
  [TMg,Fgr]=sub_read_Greenland_v3; % Greenland total FWF, km3/mo
  % Calculate annual anomalies
  % and cumulative up to date
  DVg = datevec(TMg);
  ngr = length(TMg);
  nyr =ngr/12;
  dmm = reshape(Fgr,[12,nyr]);

  Fyr = sum(dmm); % km3/yr
  Ygr = [DVg(1,1):DVg(end,1)];
  ii=find(Ygr==1990);
  Fmn = mean(Fyr(1:ii));
  ism = find(Ygr==dv0(1,1));
  cFWF = cumsum(Fyr-Fmn);
  fwf0 = cFWF(ism); % km3
  save(frv,'cFWF','Ygr');
else
  fprintf('f_griv %i, Loading %s\n',f_griv,frv);
  load(frv);
%  ii=find(Ygr==1990);
  ism = find(Ygr==dv0(1,1));
  fwf0 = cFWF(ism); % km3
end  



% Estimate volume of Greenland surplus FW in grid cell
% = ratio of FW tracer*total FWF (integrated over time 1990-date requested)
Vfw = fwf0*rr*1e9; % m3 
Hfw = Vfw./Acell; % m of FW
if IntSrf~=1
  dzRef = 50; %  rescale to m of FW in 50 m of water, for comparison
  Hfw = Hfw/hZ*dzRef;  % rescale to m of FW in 50 m of water, for comparison
else
  dzRef=hZ;
end

  
% Estiamte dlt(S) due to Greenland FWFlux anomaly
Vol = Acell.*hZ;
dS=-(Sav-(Vol-Vfw).*Sav./Vol); % freshening -> should be negative
dS(abs(dS)<1e-23)=nan;

xlim1 = 250;
xlim2 = 1210;
ylim1 = 150;
ylim2 = 1150;


% Plot fraction of the tracer 
cff=1e7;
%R=rr*cff;
R=rr;
R(R<=1e-20)=nan;
R=log(R);

btx = 'dSregions_map_v2.m';

% ------------
% Fraction of tracer
% ------------
switch(ilv)
 case(1)
  c1=-18;
  c2=-11;
 case(2)
  c1=-18;
  c2=-11;
 case(3)
  c1=-18;
  c2=-11;
 case(4)
  c1=-18;
  c2=-11;
end
nf = 1;
%stl = sprintf('Tracer Fraction *%5.3d, %i- %im, %i/%2.2i/%2.2i',...
%	      cff,abs(zz1),abs(zz2),dv0(1:3));
stl = sprintf('Ln(Tracer Fraction), %i- %im, %i/%2.2i/%2.2i',...
	      abs(zz1),abs(zz2),dv0(1:3));
if pfg1==1
  sub_plot_tracers(R,nf,HH,xlim1,xlim2,...
	 ylim1,ylim2,LON,LAT,stl,'c1',c1,'c2',c2,'cmp',5);
  bottom_text(btx,'pwd',1);
end

% ------------
% FW content, m
% ------------
switch(ilv)
 case(1)
  c1=0;
  c2=0.1;
 case(2)
  c1=0;
  c2=0.1;
 case(3)
  c1=0;
  c2=0.1;
 case(4)
  c1=0;
  c2=0.1;
end
if IntSrf==1
  c1=0;
  c2=0.5;
end

nf = 2;
stl = sprintf('GreenlAnom. dlt(FWCont), m/%im, %i- %im, %i/%2.2i/%2.2i',...
	      dzRef,abs(zz1),abs(zz2),dv0(1:3));

if pfg2==1
  sub_plot_tracers(Hfw,nf,HH,xlim1,xlim2,...
       ylim1,ylim2,LON,LAT,stl,'c1',c1,'c2',c2,'cmp',6);
  bottom_text(btx,'pwd',1);
end
% ------------
% dS change
% ------------
switch(ilv)
 case(1)
  c1=-0.1;
  c2=0;
 case(2)
  c1=-0.1;
  c2= 0;
 case(3)
  c1=-0.1;
  c2=0;
 case(4)
  c1=-0.1;
  c2=0;
end
if IntSrf==1
  c1=-0.05;
  c2=0;
end

nf = 3;
stl = sprintf('GreenlAnom. dS %i- %i, %i/%2.2i/%2.2i',...
	      abs(zz1),abs(zz2),dv0(1:3));

if pfg3==1
  sub_plot_tracers(dS,nf,HH,xlim1,xlim2,...
	 ylim1,ylim2,LON,LAT,stl,'c1',c1,'c2',c2,'cmp',3);
  bottom_text(btx,'pwd',1);
end

    
