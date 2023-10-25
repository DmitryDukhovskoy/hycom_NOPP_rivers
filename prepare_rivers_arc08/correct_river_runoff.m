% Fix some bugs in river runoff 
% (see: check_river_runoff.m)
% 1) Place Yenisey correctly
% 2) add Pyasina R.
% 3) Correct Pechora and Olenek
% Units in hycom river*.[ab] is  m/s
%c --- initialize input of river (precip bogas) forcing field
%c --- units of rivers are m/s    (positive into ocean)
%c --- rivers is always on the p grid.
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
startup;

close all
clear

PTH.data='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/force/rivers/';
PTH.topo='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
%flt=[pth,'regional.grid.b'];

ntopo=9;

%YY=2004;
flriva=sprintf('%srivers_09.a',PTH.data);
flrivb=sprintf('%srivers_09.b',PTH.data);
%flriva=sprintf('%srivers_09_Greenland_%i.a',PTH.data,YY);
%flrivb=sprintf('%srivers_09_Greenland_%i.b',PTH.data,YY);
fltopo=sprintf('%sdepth_ARCc0.08_%2.2i.nc',PTH.topo,ntopo);

flriv_new_a=sprintf('%srivers_%2.2i_corrected.a',PTH.data,ntopo); 
flriv_new_b=sprintf('%srivers_%2.2i_corrected.b',PTH.data,ntopo);

% Get topo and grid:
HH   = nc_varget(fltopo,'Bathymetry');
alat = nc_varget(fltopo,'Latitude');
elon = nc_varget(fltopo,'Longitude');
LAT  = alat;
LON  = elon;
[mm,nn]= size(HH);
[m,n]= size(HH);
IDM=n;
JDM=m;
IJDM=IDM*JDM;
npad=4096-mod(IJDM,4096);
toto=ones(npad,1);

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

% Specify rivers
% in terms of HYCOM ARCc grid;
% Conversion runoff: 0.0315* m3/s -> km3/yr
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
RV(ir).I=[147,1530];
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
%RV(ir).I=[1450,1595];
%RV(ir).J=[1210,1255];
RV(ir).I=[1450,1460];
RV(ir).J=[1210,1255];
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

% Open files with Rivers for IO:
faold = fopen(flriva,'r','ieee-be');
fbold = fopen(flrivb,'r');
fanew = fopen(flriv_new_a,'w');
fbnew = fopen(flriv_new_b,'wt');

% Write heading:
for nl=1:5
  aa=fgetl(fbold);
  disp(aa);

  if nl==4
    aa=sprintf('Corrected Yenisei, Pechora, Pyasina');
  end

  fprintf(fbnew,[aa,'\n']);
end


for k=1:12
% Read *.b:  
  aa=fgetl(fbold);
  disp(aa);

% Read HYCOM runoff:
  A=fread(faold,IJDM,'float32'); % read 2D field
  dm1=fread(faold,npad,'float32');  % Padding = size(toto)
  if size(dm1) ~= size(toto)
    error('Padding in HYCOM file ???');
  end
  toto=dm1;
  clear dm1

%  I=find(A>1e10);  % there are no land marks in rivers
%  A(I)=NaN;
  A=reshape(A,IDM,JDM)';

% Correct Pechora:
  nm0='Pechora';
  ir0=0;
  for ir=1:NR
    nm=RV(ir).Name;
    if strncmp(nm,nm0,5)
      ir0=ir;
      break;
    end
  end
  ir=ir0;
  rclim=RV(ir).Clim_m3s;
  i1=RV(ir).I(1);
  i2=RV(ir).I(2);
  j1=RV(ir).J(1);
  j2=RV(ir).J(2);
  rflx_hycom=A(j1:j2,i1:i2); % m/s
% get weight;
  wght=rflx_hycom./nansum(nansum(rflx_hycom));
  acll=ACell(j1:j2,i1:i2);
  
% Climatology  
  rcl=RV(ir).Clim_m3s(k);
  rflx=rcl*wght;
  rflx_ms=rflx*1./(acll); %m3/s-> m/s
  
  A(j1:j2,i1:i2)=rflx_ms;

  
% Correct Olenek:
  nm0='Olenek';
  ir0=0;
  for ir=1:NR
    nm=RV(ir).Name;
    if strncmp(nm,nm0,5)
      ir0=ir;
      break;
    end
  end
  ir=ir0;
  rclim=RV(ir).Clim_m3s;
  i1=RV(ir).I(1);
  i2=RV(ir).I(2);
  j1=RV(ir).J(1);
  j2=RV(ir).J(2);
  rflx_hycom=A(j1:j2,i1:i2); % m/s
% get weight;
  wght=rflx_hycom./nansum(nansum(rflx_hycom));
  acll=ACell(j1:j2,i1:i2);
  
% Climatology  
  rcl=RV(ir).Clim_m3s(k);
  rflx=rcl*wght;
  rflx_ms=rflx*1./(acll); %m3/s-> m/s
  
  A(j1:j2,i1:i2)=rflx_ms;

  
% Relocate Yenisey:
% distribute it evenly along the Yenisey channel
  nm0='Yenisey';
  ir0=0;
  for ir=1:NR
    nm=RV(ir).Name;
    if strncmp(nm,nm0,5)
      ir0=ir;
      break;
    end
  end
  ir=ir0;
  rclim=RV(ir).Clim_m3s;
  i1=RV(ir).I(1);
  i2=RV(ir).I(2);
  j1=RV(ir).J(1);
  j2=RV(ir).J(2);
  
  hh=HH(j1:j2,i1:i2);
  I=find(hh<0); % river points
  acll=ACell(j1:j2,i1:i2);
  rcl=RV(ir).Clim_m3s(k);
  wght=hh*0;
  wght(I)=1/length(I);
  rflx=rcl*wght;
  rflx_ms=rflx*1./(acll); %m3/s-> m/s
  
  A(j1:j2,i1:i2)=rflx_ms;
  
%
% Pyasina
  nm0='Pyasina';
  ir0=0;
  for ir=1:NR
    nm=RV(ir).Name;
    if strncmp(nm,nm0,5)
      ir0=ir;
      break;
    end
  end
  ir=ir0;
  rclim=RV(ir).Clim_m3s;
  i1=RV(ir).I(1);
  i2=RV(ir).I(2);
  j1=RV(ir).J(1);
  j2=RV(ir).J(2);
  rflx_hycom=A(j1:j2,i1:i2); % m/s
% get weight;
  wght=rflx_hycom./nansum(nansum(rflx_hycom));
  acll=ACell(j1:j2,i1:i2);
  
% Climatology  
  rcl=RV(ir).Clim_m3s(k);
  rflx=rcl*wght;
  rflx_ms=rflx*1./(acll); %m3/s-> m/s
  
  A(j1:j2,i1:i2)=rflx_ms;

  fprintf('Writing HYCOM files %s\n',flriv_new_a);  
  fprintf('Writing HYCOM files %s\n\n',flriv_new_b);  

  dmm=A';
  dmm=reshape(dmm,IJDM,1);
  minh=min(dmm);
  maxh=max(dmm);
  is=strfind(aa,'=');
  ch=aa(1:is);
  anew=sprintf('%s%3i    %10.7E   %10.7E',ch,k,minh,maxh);
  fprintf(fbnew,[anew,'\n']);
  fwrite(fanew,dmm,'float32','ieee-be');
  fwrite(fanew,toto,'float32','ieee-be');

end;  % for k - months

fclose(fanew);
fclose(fbnew);
fclose(faold);
fclose(fbold);











