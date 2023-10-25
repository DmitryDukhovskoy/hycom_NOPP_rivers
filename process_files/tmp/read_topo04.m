% Read binary topo file
%
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close

PTH.topo='/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';
T = 17;            % topo #
pthT = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';
%pthout=sprintf('/Net/bimonthly/ddmitry/%s/%s/topo_nc/',RG,EX);
pthout=pthT;

% Topo/grid files, input:
%flda=sprintf('%sdepth_ARCc0.04_%2.2i.a',pthT,T);
flda=sprintf('%sdepth_ARCc0.04_%2.2iDD.a',pthT,T);
%flda=sprintf('%sdepth_ARCc0.04_11D.a',pthT);
flga=sprintf('%sregional.grid.a',pthT);
flgb=sprintf('%sregional.grid.b',pthT);

fprintf('Reading %s\n',flda);
fprintf('Reading %s\n',flga);

% Read in topo from *a file:
f1=fopen(flgb);  % read I,J from regional.grid.b
aa=fgetl(f1);
dmm=aa(2:8);
IDM=str2num(dmm);
aa=fgetl(f1);
dmm=aa(2:8);
JDM=str2num(dmm);
IJDM=IDM*JDM;
fclose(f1);

npad=4096-mod(IJDM,4096);

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % read lon/lat from regional grid file
grid_fid=fopen(flga,'r');

[plon,count]=fread(grid_fid,IJDM,'float32','ieee-be');
stat=fseek(grid_fid,4*(npad+IJDM),-1);
if stat<0
  error('Reading grid file ...');
end
[plat,count]=fread(grid_fid,IJDM,'float32','ieee-be');

disp('Reading lat/lon  ...')
plon=(reshape(plon,IDM,JDM))';
plat=(reshape(plat,IDM,JDM))';

fclose(grid_fid);

% --------------------------
I=find(plon>180);
plon(I)=plon(I)-360;
I=find(plon<-180);
plon(I)=plon(I)+360;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read bathymetry from regional.depth.a

depth_fid=fopen(flda,'r');

%fseek(depth_fid,6*4*(npad+IJDM),-1) % this is confusing, should not be needed
[h,count]=fread(depth_fid,IJDM,'float32','ieee-be');
%y=find(h>1e10);
%h(y)=nan;
HH=reshape(h,IDM,JDM)';
%clear h;
fclose(depth_fid);

HH(HH>1e10)=nan;
HH=-HH;
HH(HH>-0.1)=0;
HH(isnan(HH))=100;




