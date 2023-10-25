% Check iso.sigma file - varyaing isopycnal layers
% ARCc0.08 - 32 or 41 layers
% for ARCc0.04 T17DD (corrected T17) layers 41 
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

%kold=41;
knew=41;
%IDM=1600; % ARCc0.08 (old) grid
%JDM=2520; % ARCc0.08
%IJDM=IDM*JDM;
%npad=4096-mod(IJDM,4096);
%toto=ones(npad,1);

% Read in sigma layers:
%fina  = sprintf('%siso_sigma_41lr.a',pthin);
%finb  = sprintf('%siso_sigma_41lr.b',pthin);
fouta = sprintf('%siso_sigma_T%s_L%2.2i.a',pthout,TV,knew);
foutb = sprintf('%siso_sigma_T%s_L%2.2i.b',pthout,TV,knew);

fprintf('Opening %s\n',fouta);
fprintf('Opening %s\n',foutb);
fida = fopen(fouta,'r','ieee-be');


% Get new topo and grid:
fltopo_new=sprintf('%sdepth_ARCc0.04_%s.nc',pthtopo,TV);
HH   = nc_varget(fltopo_new,'Bathymetry');
alat = nc_varget(fltopo_new,'Latitude');
elon = nc_varget(fltopo_new,'Longitude');
LAT  = alat;
LON  = elon;
[mm,nn]= size(HH);
IDM  = nn;
JDM  = mm;
IJDM = IDM*JDM;
npad = 4096-mod(IJDM,4096);
toto = ones(npad,1);


% Target dneisities in the 41-layer version
pthb='/home/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers/';
flblk=sprintf('%sblkdat.input_ARCc0.08_41lev',pthb);
TDENS=read_targ_dens_blkdat(flblk);


for ll=1:knew
  A=fread(fida,IJDM,'float32'); % read 2D field
  dm1=fread(fida,npad,'float32');  % Padding = size(toto)
  minh=min(A);
  maxh=max(A);

  fprintf('Reading old iso_sigma, k=%i, minh/maxh: %9.4f %9.4f\n',ll,minh,maxh);

  f_chck=10;
  if ll==f_chck
    AA=reshape(A,IDM,JDM)';
    figure(10); clf;
    nl=ll;
    axes('Position',[0.1 0.1 0.8 0.8]);
    hold on;
    pcolor(AA); shading flat
    contour(HH,[0 0],'k');
    title('New iso_sigma');
    caxis([0 10]);
    colorbar

  end
%keyboard

end

fclose(fida);


