% CHeck river runoff in the 
% hycom river_*.[ab] files
% Compare with AOMIP climatology
% Units:  m/s
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
startup;

close all
clear

PTH.data='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/force/rivers/';
PTH.topo='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
PTH.dataout = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/force/riversNCAR/';
%flt=[pth,'regional.grid.b'];

TT=11; % topography

YY=2016;
%YY=2004;
%flriva=sprintf('%srivers_09.a',PTH.data);
%flrivb=sprintf('%srivers_09.b',PTH.data);
%flriva=sprintf('%srivers_09_corrected.a',PTH.data);
%flrivb=sprintf('%srivers_09_corrected.b',PTH.data);
flriva=sprintf('%srivers_%2.2i_NCAR_Gr_%i.a',PTH.dataout,TT,YY);
flrivb=sprintf('%srivers_%2.2i_NCAR_Gr_%i.b',PTH.dataout,TT,YY);
fltopo=sprintf('%sdepth_ARCc0.08_%2.2i.nc',PTH.topo,TT);

%fltxt=[pth1,fvrs,'.b'];
%flbth=[pth1,fvrs,'.a'];
%flrv=[pth,'rivers_03.a'];

% Get topo and grid:
HH   = nc_varget(fltopo,'Bathymetry');
alat = nc_varget(fltopo,'Latitude');
elon = nc_varget(fltopo,'Longitude');
LAT  = alat;
LON  = elon;
[mm,nn]= size(HH);
[m,n]= size(HH);

DX=zeros(mm,nn);
DY=zeros(mm,nn);
for i=1:nn-1
  dx=distance_spheric_coord(LAT(:,i),LON(:,i),LAT(:,i+1),LON(:,i+1));
  DX(:,i)=dx;
end
DX(:,nn)=dx;
for j=1:mm-1
  dy=distance_spheric_coord(LAT(j,:),LON(j,:),LAT(j+1,:),LON(j+1,:));
  DY(j,:)=dy;
end
DY(mm,:)=dy;

ACell=DX.*DY;



% Reading rivers & bathymetry:
IDM=n;
JDM=m;
IJDM=IDM*JDM;
npad=4096-mod(IJDM,4096);

% Specify rivers
% in terms of HYCOM ARCc grid;
% Conversion runoff: 0.0315* annual mean(m3/s) -> km3/yr
% Annual gauged vol. flux is 2456 km3/yr (AWI climatology)
% AOMIP river runoff data 
ir=0;

ir=ir+1;
RV(ir).Name='Mackenzie';
RV(ir).I=[400,450];
RV(ir).J=[1600,1700];
RV(ir).Clim_m3s=[3814.35, 3605.9, 3349.0, 3384.45,...
		 13218.8, 21413.0, 17854.7, 13984.2, ...
		 11268.7, 9038.15, 4756.0, 3595.0];

ir=ir+1;
RV(ir).Name='S.Dvina';
RV(ir).I=[1470,1530];
RV(ir).J=[650,800];
RV(ir).Clim_m3s=[1033.45, 826.6, 724.9, 2415.1, ...
		 13839.4, 7029.7, 2943.8, 2149.4,...
		 2320.4, 2912.4, 2363.0, 1400.7];

ir=ir+1;
RV(ir).Name='Pechora';
RV(ir).I=[1490,1530];
RV(ir).J=[900,1000];
RV(ir).Clim_m3s=[959.4, 773.8, 695.4, 950.2,...
		 15502.6, 17126.4, 5534.2, 3227.8,...
		 3917.0, 4197.6, 1894.4, 1277.4];

ir=ir+1;
RV(ir).Name='Ob';
RV(ir).I=[1520,1590];
RV(ir).J=[1110,1200];
RV(ir).Clim_m3s=[4986.7, 4120.2, 3635.8, 3698.0,...
		 15122.7, 36715.9, 31694.8, 23122.7,...
		 14747.4, 11000.7, 6695.4, 5734.2];

ir=ir+1;
RV(ir).Name='Yenisey';
RV(ir).I=[1450,1600];
RV(ir).J=[1200,1260];
RV(ir).Clim_m3s=[6038.6, 6022.9, 5983.9, 6001.3, ...
		 27533.5, 77386.6, 26586.8, 17485.4,...
		 16896.1, 13969.2, 6855.9, 5839.8];

ir=ir+1;
RV(ir).Name='Pyasina';
RV(ir).I=[1400,1500];
RV(ir).J=[1300,1350];
RV(ir).Clim_m3s=[501.1, 501.1, 501.1, 501.1,...
		 501.1, 7516.7, 10022.3, 2818.8,...
		 3758.4, 4071.6, 1252.8, 689.0];

ir=ir+1;
RV(ir).Name='Olenek';
RV(ir).I=[1300,1380];
RV(ir).J=[1550,1600];
RV(ir).Clim_m3s=[7.00000, 3.20000, 2.10000, 1.65000,...
		 306.850, 7160.00, 2203.50, 846.100,...
		 1050.45, 308.900, 78.7500, 25.1000];

ir=ir+1;
RV(ir).Name='Lena';
RV(ir).I=[1150,1300];
RV(ir).J=[1600,1750];
RV(ir).Clim_m3s=[2783.04, 2136.78, 1651.78, 1350.20,...
		 6235.90, 73917.4, 39683.4, 27340.0,...
		 24126.0, 13771.7, 3502.06, 2927.86];

ir=ir+1;
RV(ir).Name='Kolyma';
RV(ir).I=[930,1000];
RV(ir).J=[1840,1870];
RV(ir).Clim_m3s=[260.0, 191.765, 190.706, 162.647,...
		 1739.24, 13996.6, 7527.12, 6011.70,...
		 4804.71, 1646.18, 436.059, 340.647];

NR=ir;

% Rivers:
riv_fid=fopen(flriva,'r','ieee-be');

for k=1:12  % time - 12mo
  fprintf('check_river_runoff:    Reading month %i\n',k);

  A  = fread(riv_fid,IJDM,'float32'); % read 2D field
  dm1= fread(riv_fid,npad,'float32');  % Padding = size(toto)

  I=find(A>1e10);
  A(I)=NaN;
  A=reshape(A,IDM,JDM)';
  A(A==0)=nan;
  
% Check max/min:
%keyboard
  fprintf('*a: min=%8.5d max=%8.5d\n',min(min(A)),max(max(A)));

% Convert m/s -> m3/s  
  RFlx=A.*ACell;

% Get Runoff by rivers
  for ir=1:NR
    i1=RV(ir).I(1);
    i2=RV(ir).I(2);
    j1=RV(ir).J(1);
    j2=RV(ir).J(2);
    dmm=RFlx(j1:j2,i1:i2);
    Rtot=nansum(nansum(dmm));
    RV(ir).R_hycom_m3s(k)=Rtot;
  end

  JT=find(~isnan(RFlx));
  IL=find(HH(JT)>=0);
  if ~isempty(IL),
    error('check_river_runoff.m: River source is ONLAND ');
  end
  
end

% Plot
Ftot=zeros(12,1);
Fclim=zeros(12,1);
figure(1);
for ir=1:NR
  nm=RV(ir).Name;
  rcl=RV(ir).Clim_m3s;
  rhycom=RV(ir).R_hycom_m3s;
  
  subplot(3,3,ir);
  plot(rcl,'r','linewidth',2);
  hold on;
  plot(rhycom,'b','linewidth',2);
  stt=sprintf('%s, m3/s, hycom(b)',nm);
  title(stt);
  set(gca,'tickdir','out',...
	  'xlim',[1 12]);
  Ftot=Ftot+rhycom';
  Fclim=Fclim+rcl';
end;

% Annual mean: km3/yr
Fyr=mean(Ftot)*0.0315;
Fyr_cl=mean(Fclim)*0.0315; % some rivers are missing here
fprintf('Annual runoff, km3/yr in HYCOM: %8.3f, Clim=2456\n',Fyr);

txtbtm='hycom_arc08/prepare_rivers/check_river_runoff.m';
bottom_text(txtbtm);







