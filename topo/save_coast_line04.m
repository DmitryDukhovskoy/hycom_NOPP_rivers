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
pthout = pthT;
fcst   = sprintf('%scoast_line_arc04.dat',pthout); 

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


% Specify region for coast:
f_cst = 'Laptev';
switch(f_cst),
 case('Laptev');
  IC = [1866,2734,2734,1866];
  JC = [3593,3593,2675,2675];
end
i1 = min(IC);
j1 = min(JC);
i2 = max(IC);
j2 = max(JC);
cH = HH(j1:j2,i1:i2);
cX = plon(j1:j2,i1:i2);
cY = plat(j1:j2,i1:i2);

figure(2); clf;
[cc1,cc2] = contour(cX,cY,cH,[0 0],'b');
hold on;


np=size(cc1,2);
n0=1000;    % min # of points in the LC contour
ic=0;
kk=1;
clear COAST
while kk<np
  nrd=cc1(2,kk);
  iiS=kk+1;
  iiE=kk+nrd;
  xx=cc1(1,iiS:iiE)';
  yy=cc1(2,iiS:iiE)';

  chck=0;
  if nrd>n0
    ic=ic+1;
    COAST(ic).xx=xx;
    COAST(ic).yy=yy;
  end
  kk=iiE+1;
end; % while kk<np  

% Pick the longest
if ic>1
  amx = 0;
  for k=1:ic
    xx=COAST(k).xx;
    lng(k,1) = length(xx);
  end
  
  imx = find(lng == max(lng));
else
  imx = 1;
end

XX = COAST(imx).xx;
YY = COAST(imx).yy;
plot(XX,YY,'r.');

fout = sprintf('%scoast_%s.dat',pthout,f_cst);
fprintf('Writing data to %s\n',fout);
fid  = fopen(fout,'w');

nn = length(XX);
for k=1:nn
  aa = sprintf('%11.7f %11.7f',XX(k),YY(k));
  fprintf(fid,[aa,'\n']);
end
fclose(fid);














