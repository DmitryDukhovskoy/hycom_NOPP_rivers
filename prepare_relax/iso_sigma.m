% Create iso.sigma file - varyaing isopycnal layers
% for ARCc T11 layers 41 
%
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close

ptharc  = '/Net/kronos/ddmitry/hycom/ARCc0.08/tmp_files/';
%pthglb  = '/Net/kronos/ddmitry/hycom/GLBb0.08/expt_69.1/restart/';
pthin  = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/relax/010/';
pthout = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/relax/110/41layers_T11/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';

ntopo1=09;
ntopo2=11;

% Read in sigma layers:
fina  = sprintf('%siso_sigma.a',pthin);
finb  = sprintf('%siso_sigma.b',pthin);
fouta = sprintf('%siso_sigma_41lr.a',pthout);
foutb = sprintf('%siso_sigma_41lr.b',pthout);

fida = fopen(fina,'r','ieee-be');

kold=32;
knew=41;
IDM=1600;
JDM=2520;
IJDM=IDM*JDM;
npad=4096-mod(IJDM,4096);
toto=ones(npad,1);


% Get new topo and grid:
fltopo_new=sprintf('%sdepth_ARCc0.08_%2.2i.nc',pthtopo,ntopo2);
HH   = nc_varget(fltopo_new,'Bathymetry');
alat = nc_varget(fltopo_new,'Latitude');
elon = nc_varget(fltopo_new,'Longitude');
LAT  = alat;
LON  = elon;
[mm,nn]= size(HH);


fltopo_old=sprintf('%sdepth_ARCc0.08_%2.2i.nc',pthtopo,ntopo1);
HHo  = nc_varget(fltopo_old,'Bathymetry');

% Target dneisities in the 41-layer version
pthb='/home/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers/';
flblk=sprintf('%sblkdat.input_ARCc0.08_41lev',pthb);
TDENS=read_targ_dens_blkdat(flblk);


for ll=1:kold
  fprintf('Reading old iso_sigma, k=%i\n',ll);
  A=fread(fida,IJDM,'float32'); % read 2D field
  dm1=fread(fida,npad,'float32');  % Padding = size(toto)

  A=reshape(A,IDM,JDM)';

  P1(ll,:,:)=A;
end

fclose(fida);


% Baltic Sea uses different sigma-layers: z-layers
% create dz array
dz=[1:knew]'*0.1;
A=squeeze(P1(1,:,:));
I=find(A<10); % where hybrid/isopycnic layers changed to fixed z
Pnew=zeros(knew,JDM,IDM);
Pnew=repmat(TDENS,[1 JDM IDM]);

for k=1:knew
  Pnew(k,I)=dz(k);
end

f_chck=0;
if f_chck>0
  figure(10); clf;
  nl=1;
  A=squeeze(P1(nl,:,:));
  axes('Position',[0.05 0.51 0.6 0.45]);
  contour(HH,[0 0],'k');
  hold on;
  pcolor(A); shading flat
  title('Old iso_sigma');
  caxis([0 10]);

  A=squeeze(Pnew(nl,:,:));
  axes('Position',[0.05 0.03 0.6 0.45]);
  contour(HH,[0 0],'k');
  hold on;
  pcolor(A); shading flat
  title('New iso_sigma');
  caxis([0 10]);
  
end

%
% Write output
fidaOUT = fopen(fouta,'w');
fidbOUT = fopen(foutb,'wt');

fprintf('Writing HYCOM files %s\n',fouta);  
fprintf('Writing HYCOM files %s\n\n',foutb);  

cline='target_density: layer,range =';
for kk=1:knew
  A = squeeze(Pnew(kk,:,:));
  dmm=A';
  dmm=reshape(dmm,IJDM,1);
  minh=min(dmm);
  maxh=max(dmm);
  aa=sprintf('%s %2.2i %9.4f %9.4f',cline,kk,minh,maxh);
  fprintf('%s\n',aa);
  fprintf(fidbOUT,[aa,'\n']);
  fwrite(fidaOUT,dmm,'float32','ieee-be');
  fwrite(fidaOUT,toto,'float32','ieee-be');
end

fclose(fidaOUT);
fclose(fidbOUT);

