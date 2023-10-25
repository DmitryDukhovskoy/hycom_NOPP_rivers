% Check  veldf4
% for ARCc0.04 T17DD (corrected T17) layers 41 
% Alan: 
% Dmitry,
%
%veldf4 has already been increased to 0.02  in this region:
%
%narwhal04 3611> head *df4.b
%==> thkdf4.b <==
%thkdf4: range =    0.0100    0.0100
%
%==> veldf4.b <==
%veldf4: range =    0.0100    0.0200
%
%narwhal04 3613> hycom_range veldf4.a 3200 5040 1976 629
%min, max =   9.99999419E-03  1.99999996E-02    (1976, 629) =   1.99999996E-02
%
% You could try increasing veldf4 in the same patch to 0.04 (say) .  Note that there is a CFL limit on veldf4, so don't increase it too much.
%
% In both GLBc0.04 and ATLc0.02 we have increased thkdf4 as high as 0.07, and in ATLc0.02 this might be in the region where you are blowing up.
%
% If you do increase a df4 field, it is best to repeat the previous month to give it time to work.
%
% Finally, it is best to set batrop as large as possible (consistent with baclin).  
% You are using 2.5=baclin/20, but baclin/14=3.571428 or baclin/16=3.125 should work. 
% If batrop is too large the model will blow up almost immediately, so it is easy to test for its maximum value.

addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close


pthdf = '/Net/mars/ddmitry/hycom/ARCc0.04/relax/023/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';

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
%fouta = sprintf('%sveldf4.a',pthdf);
%foutb = sprintf('%sveldf4.b',pthdf);
fouta = sprintf('%sveldf4_0406_T17DD.a',pthdf);
foutb = sprintf('%sveldf4_0406_T17DD.b',pthdf);

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


  A=fread(fida,IJDM,'float32'); % read 2D field
  dm1=fread(fida,npad,'float32');  % Padding = size(toto)
  minh=min(A);
  maxh=max(A);

  fprintf('Reading %s, minh/maxh: %9.4f %9.4f\n',fouta,minh,maxh);

  f_chck=1;
  if f_chck>0
    AA=reshape(A,IDM,JDM)';
    figure(10); clf;
    axes('Position',[0.1 0.1 0.8 0.8]);
    hold on;
    pcolor(AA); shading flat
    title(sprintf('%s',fouta));
    caxis([0 maxh]);
    contour(HH,[0 0],'k');
    colorbar
    
  end
%keyboard


fclose(fida);


