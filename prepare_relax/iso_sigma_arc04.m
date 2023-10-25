% Create iso.sigma file - varyaing isopycnal layers
% for ARCc0.04 T17DD (corrected T17) layers 41 
% Need to interpolate from ARCc0.08 into 0.04
% For low-saline regions (Baltic Sea and usually Black Sea - if exists)
% isop. layers are set to fixed Z
% need to specify dz in these regions
%
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close


pthin   = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/relax/110/41layers_T11/';
pthout  = '/Net/mars/ddmitry/hycom/ARCc0.04/relax/010/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';

ntopo1=11;
ntopo2=17;
TV = sprintf('%2.2iDD',ntopo2);

kold=41;
knew=41;
IDM=1600; % ARCc0.08 (old) grid
JDM=2520; % ARCc0.08
IJDM=IDM*JDM;
npad=4096-mod(IJDM,4096);
%toto=ones(npad,1);



% Read in sigma layers:
fina  = sprintf('%siso_sigma_41lr.a',pthin);
finb  = sprintf('%siso_sigma_41lr.b',pthin);
fouta = sprintf('%siso_sigma_T%s_L%2.2i.a',pthout,TV,knew);
foutb = sprintf('%siso_sigma_T%s_L%2.2i.b',pthout,TV,knew);

fida = fopen(fina,'r','ieee-be');


% Get new topo and grid:
fltopo_new=sprintf('%sdepth_ARCc0.04_%s.nc',pthtopo,TV);
HH   = nc_varget(fltopo_new,'Bathymetry');
alat = nc_varget(fltopo_new,'Latitude');
elon = nc_varget(fltopo_new,'Longitude');
LAT  = alat;
LON  = elon;
[mm,nn]= size(HH);
IDMn  = nn;
JDMn  = mm;
IJDMn = IDMn*JDMn;
npadN = 4096-mod(IJDMn,4096);
totoN = ones(npadN,1);

% Target dneisities in the 41-layer version
pthb='/home/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers/';
flblk=sprintf('%sblkdat.input_ARCc0.08_41lev',pthb);
TDENS=read_targ_dens_blkdat(flblk);


for ll=1:kold
  A=fread(fida,IJDM,'float32'); % read 2D field
  dm1=fread(fida,npad,'float32');  % Padding = size(toto)
  minh=min(A);
  maxh=max(A);

  fprintf('Reading old iso_sigma, k=%i, minh/maxh: %9.4f %9.4f\n',ll,minh,maxh);

  A=reshape(A,IDM,JDM)';

  P1(ll,:,:)=A;
end

fclose(fida);

%keyboard

% Baltic Sea uses different sigma-layers: z-layers
% create dz array
dz=[1:knew]'*0.1;
A=HH*0;
dmm=HH(540:1250,2561:2911);
A(540:1250,2561:2911)=dmm; % where hybrid/isopycnic layers changed to fixed z
I=find(A<0); 
Pnew=zeros(knew,JDMn,IDMn);
Pnew=repmat(TDENS,[1 JDMn IDMn]);

for k=1:knew
  fprintf('iso_sigma, fixed z: k=%2.2i\n',k);
  Pnew(k,I)=dz(k);
end

f_chck=0;
if f_chck>0
  figure(10); clf;
  nl=1;
  A=squeeze(P1(nl,:,:));
  axes('Position',[0.05 0.51 0.6 0.45]);
%  contour(HH,[0 0],'k');
  hold on;
  pcolor(A); shading flat
  title('Old iso_sigma');
  caxis([0 10]);

  A=squeeze(Pnew(nl,:,:));
  axes('Position',[0.05 0.03 0.6 0.45]);
  hold on;
  pcolor(A); shading flat
  contour(HH,[0 0],'k');
  title('New iso_sigma');
  caxis([0 10]);
  
end

%keyboard

%
% Write output
fidaOUT = fopen(fouta,'w');
fidbOUT = fopen(foutb,'wt');

fprintf('Writing HYCOM files %s\n',fouta);  
fprintf('Writing HYCOM files %s\n\n',foutb);  

cline='target_density: layer,range =';
for kk=1:knew
  A = squeeze(Pnew(kk,:,:));
  dmm = A';
  dmm = reshape(dmm,IJDMn,1);
  minh = min(dmm);
  maxh = max(dmm);
  aa = sprintf('%s %2.2i %9.4f %9.4f',cline,kk,minh,maxh);
  fprintf('%s\n',aa);
  fprintf(fidbOUT,[aa,'\n']);
  fwrite(fidaOUT,dmm,'float32','ieee-be');
  fwrite(fidaOUT,totoN,'float32','ieee-be');
%keyboard
end

fclose(fidaOUT);
fclose(fidbOUT);
fprintf('iso.sigma ready: %s\n', fouta);
fprintf('iso.sigma ready: %s\n', foutb);
