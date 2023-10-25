% Spatial map of dltS for specified day
%
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

dnmb = datenum(2016,12,30);
dv0  = datevec(dnmb);
yr   = dv0(1);
iday = dnmb-datenum(yr,1,1)+1;

nTr = 1; 
%ilv = 1; % 0-50m
%ilv = 2; % 50-150m
%ilv = 3; % bottom-surf
ilv = 4; % 150-300m

LR(1,:) = [0,-50];
LR(2,:) = [-50,-150];
LR(3,:) = [0,-10000]; % full water column
LR(4,:) = [-150, -300];
nlrs= 4; 

zz1 = LR(ilv,1);
zz2 = LR(ilv,2);

dz=abs(zz2-zz1);

fprintf('S change, ilv=%i, %i - %i, nTr=%i, %i/%i/%i\n',ilv,zz1,zz2,nTr,dv0(1:3));

regn = 'ARCc0.08';
expt = 110;  
pthfig  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/fig_NAtl/',expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
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

pthbin  = sprintf('/nexsan/archive/ARCc0.08_110/data/%4.4i/',yr);
fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);
fld='salin';
[SS,n,m,l] = read_hycom(fina,finb,fld);
SS(SS>1e20)=nan;
[ZM,ZZ] = sub_zz_zm(fina, finb,HH,'f_btm',1);
% Find mean S for layer ilv
Zk  = HH*0;
cZk = HH*0; % count v. layers
for kk=1:l
  zz= squeeze(ZZ(kk,:,:));
  I=find(HH<zz2 & zz>=zz2 & zz<zz1);
  if ilv==3
    zzn= squeeze(ZZ(kk+1,:,:));
    I=find(~isnan(zz) & ~isnan(zzn)); % find bottom
  end
  if isempty(I); continue; end;
  if max(max(zz))<zz2, break; end;
  Zk(I)=kk;
  cZk(I)=cZk(I)+1;
end

Zk(Zk==0)=nan; % bottom layer within the depth interval
Zt=Zk-cZk+1;   % top layer
Sav = HH*0; % S averaged within the layer
lmx = max(max(Zk));
hZ  = HH*0;
smm = HH*0;
for kk=1:lmx
  I=find(Zk>=kk & Zt<=kk);
  if isempty(I), continue; end;
  dZ=abs(squeeze(ZZ(kk+1,:,:))-squeeze(ZZ(kk,:,:)));
  hZ(I)=hZ(I)+dZ(I);
%  if isnan(hZ(300,600)); keyboard; end;
  smm = smm+squeeze(SS(kk,:,:)).*dZ;
  Sav(I) = smm(I)./hZ(I); % depth-avrg S
end
Sav(Sav==0)=nan;


% Tracer fraction in grid cells
rr = sub_fraction_tracer(expt,regn,Acell,HH,dnmb,nTr,ilv,hZ);
%rr = sub_fraction_tracerMass(expt,regn,Acell,HH,dnmb,nTr,ilv,hZ);

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
keyboard
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
  c1=-16;
  c2=-9;
 case(4)
  c1=-18;
  c2=-11;
end
nf = 1;
%stl = sprintf('Tracer Fraction *%5.3d, %i- %im, %i/%2.2i/%2.2i',...
%	      cff,abs(zz1),abs(zz2),dv0(1:3));
stl = sprintf('Ln(Tracer Fraction), %i- %im, %i/%2.2i/%2.2i',...
	      abs(zz1),abs(zz2),dv0(1:3));
sub_plot_tracers(R,nf,HH,xlim1,xlim2,...
       ylim1,ylim2,LON,LAT,stl,'c1',c1,'c2',c2,'cmp',5);
btx = 'dSregions_map.m';
bottom_text(btx,'pwd',1);

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
  c2=1;
 case(4)
  c1=0;
  c2=0.1;
end
nf = 2;
stl = sprintf('GreenlAnom. dlt(FWCont), m, %i- %im, %i/%2.2i/%2.2i',...
	      abs(zz1),abs(zz2),dv0(1:3));
sub_plot_tracers(Hfw,nf,HH,xlim1,xlim2,...
       ylim1,ylim2,LON,LAT,stl,'c1',c1,'c2',c2,'cmp',6);
bottom_text(btx,'pwd',1);

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

nf = 3;
stl = sprintf('GreenlAnom. dS %i- %i, %i/%2.2i/%2.2i',...
	      abs(zz1),abs(zz2),dv0(1:3));
sub_plot_tracers(dS,nf,HH,xlim1,xlim2,...
       ylim1,ylim2,LON,LAT,stl,'c1',c1,'c2',c2,'cmp',3);
bottom_text(btx,'pwd',1);


    
