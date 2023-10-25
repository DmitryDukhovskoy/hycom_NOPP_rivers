% Read rivers
% Units:  m/s
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
startup;

close all
clear

PTH.data='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/force/rivers/';
PTH.topo='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
%flt=[pth,'regional.grid.b'];

YY=2004;
T=11;
flriva=sprintf('%srivers_%2.2i.a',PTH.data,T);
flrivb=sprintf('%srivers_%2.2i.b',PTH.data);
%flriva=sprintf('%srivers_09.a',PTH.data);
%flrivb=sprintf('%srivers_09.b',PTH.data);
%flriva=sprintf('%srivers_%2.2i_Greenland_%i.a',PTH.data,T,YY);
%flrivb=sprintf('%srivers_%2.2i_Greenland_%i.b',PTH.data,T,YY);
fltopo=sprintf('%sdepth_ARCc0.08_%2.2i.nc',PTH.topo,T);

fprintf('Opening %s\n',flriva);

%fltxt=[pth1,fvrs,'.b'];
%flbth=[pth1,fvrs,'.a'];
%flrv=[pth,'rivers_03.a'];

% Get topo and grid:
HH   = nc_varget(fltopo,'Bathymetry');
alat = nc_varget(fltopo,'Latitude');
elon = nc_varget(fltopo,'Longitude');
[m,n]= size(HH);

% Reading rivers & bathymetry:
IDM=n;
JDM=m;
IJDM=IDM*JDM;
npad=4096-mod(IJDM,4096);


% Rivers:
riv_fid=fopen(flriva,'r');
%fseek(riv_fid,6*4*(npad+IJDM),-1);
for im=1:12
  fprintf('Month %i\n',im);
%  dmm=fread(riv_fid,1,'int','ieee-be'); <-- Do not need this, direct access format
  [A,counta]=fread(riv_fid,IJDM,'float32','ieee-be');
  toto = fread(riv_fid,npad,'float32','ieee-be');
%  dmm=fread(riv_fid,1,'int','ieee-be'); <-- Do not need this, direct access format
%  I=find(A>1e10);
%  A(I)=0;
  A=reshape(A,IDM,JDM)';
  if im==1,
    A1=A;
  end
  
end;


[jj,ii]=find(A>0);


clf
contour(HH,[0 0],'k');
hold on;
contour(A1,[0 0],'g');
contour(A,[0 0],'r--'); % should match A1
keyboard

tmp=log10(A);
tmp(A==0)=NaN;
%pcolor(tmp); shading flat;
%caxis([-10 -2]);
%colorbar
%title('log10(River runoff), m/s')

%for j=1:length(jj)
%  i0=ii(j);
%  j0=jj(j);
%  x=tmp(j0,i0);
%  intr=max(find(cnt<=x));
%  if isempty(intr), intr=1; end;
%  p1=plot(i0,j0,'r.');
%  set(p1,'MarkerSize',6,'Color',jjs(intr,:));
%end;

% Check if rivers are in the ocean points:
B2=HH;
J=find(B2==0);
B2(J)=NaN;
J2=find(~isnan(B2));
B2(J2)=0;
B2(J)=1;

J3=find(tmp~=0);
tmp(J3)=1;

D=tmp+B2;  % 
 I=find(D>1);
if ~isempty(I); error('Some river points are on land'); end








