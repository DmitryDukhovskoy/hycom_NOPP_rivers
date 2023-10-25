% Plot tracer fluxes and tracer budget
% at the OB - extracted in trOB.m
% to analyze Tracer fluxes
% from the OB into the domain
%
% mat files are saved by years
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

yr1 = 1993;
yr2 = 2015;


fprintf('======    Years: %i-%i\n',yr1,yr2);

rg = 9806; 
%cBr = 14; % coefficient for Bering Str. tracer - not right
cBr = sub_BeringTrConc(0);% coefficient for Bering Str., tracer, kg/m3

regn = 'ARCc0.08';
expt = 110;  
%pthfig  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/fig_trac/',expt);
pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);


ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2

hmsk=HH;
hmsk(HH<0)=nan;

% Specify segments of the x-section
i1=193;
i2=1021;
j1=21;
j2=j1;
IJs=[ i1  j1;
      i2  j2];
Rmsk = HH;
Rmsk(j1+1:end,:)=0;
Rmsk(:,1:i1)=0;
Rmsk(:,i2+1:end)=0;

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

SEGM.I = IJs(:,1);
SEGM.J = IJs(:,2);
SEGM.Ind = INDs;
SEGM.dx  = DX(INDs);
SEGM.nx  = 0; % norm, x component
SEGM.ny  = 1; % norm, y component



FTr=[];
Mob=[];
Mtot=[];
TM=[];
for iyr = yr1:yr2
  yr=iyr;
  fmat = sprintf('%strcrFlx_AtlOB_daily_%4.4i.mat',pthmat,iyr);
  fprintf('Loading %s\n',fmat);
  load(fmat);
  
  FTr=[FTr;TRFLX.FLX_kg_s];
  TM=[TM;TRFLX.TM];
  Mob=[Mob;TRFLX.SpongeReg_Mass_ton]; % Tr Mass in Sponge Region
  Mtot=[Mtot;TRFLX.Overall_Mass_ton]; % Tr Mass in whole domain
end

figure(1); clf;
hold on
contour(HH,[0 0],'k');
contour(HH,[-5000:1000:-100],'Color',[0.7 0.7 0.7]);
plot(IIs,JJs,'b.-');
axis('equal');
set(gca,'xlim',[180 1200],...
	'ylim',[0 300]);



% Plot Total transport towards Atl OB;
DV = datevec(TM);
nrc = length(TM);
clear YR
for ik=1:nrc
  yr=DV(ik,1);
  nyd=365;
  if mod(yr,4)==0
    nyd=366;
  end
  YR(ik)=(TM(ik)-datenum(yr,1,1))/nyd+yr;
end

figure(2); clf
axes('position',[0.08 0.7 0.86 0.25]);
plot(YR,FTr);
set(gca,'tickdir','out',...
	'xlim',[YR(1) YR(end)],...
	'xtick',[1990:2:YR(end)],...
	'xgrid','on',...
	'ygrid','on');
title('Tracer transp (ton/s) near Atl.OB, + northward');   

% Time-integrated tracer mass flux, kg
dt=diff(TM);
dt=[dt(1);dt];
cM = cumsum(FTr.*dt*3600);  % kg
axes('position',[0.08 0.4 0.86 0.25]);
plot(YR,cM);
set(gca,'tickdir','out',...
	'xlim',[YR(1) YR(end)],...
	'xtick',[1990:2:YR(end)],...
	'xgrid','on',...
	'ygrid','on');
title('Time-integrated TrMass gain across AtlOB, kg');

%figure(3);
% Time integrated mass flux
% as a fraction to the total
% Tracer mass in the domain
axes('position',[0.08 0.08 0.86 0.25]);
plot(YR,cM./(Mtot*1e3)); % kg/kg
set(gca,'tickdir','out',...
	'xlim',[YR(1) YR(end)],...
	'xtick',[1990:2:YR(end)],...
	'xgrid','on',...
	'ygrid','on');
title('Fraction of TrMass gain across AtlOB wrt total TrMass Domain');

btx = 'plot_trOB.m';
bottom_text(btx,'pwd',1);

    

