% Interpolate NRL offlux into my T11 topography
% T11 is similar to the NRL
% no need to interpolate the offlux file
%
% However, I need to create 0 offlux file 
% to reduce sea ice melt in the Arctic
% Prepare 0-offset
%
% From Alan's email: 
% The flux offset field corrects the atmospheric fluxes in HYCOM 
% to get a better long term mean SST (and long term average ice concentration).
%  In HYCOM this is a separate input field, but in CICE we add the 
% same field to longwave downward flux as a preprocessing step.  
% Setting it to zero is OK, but you will have larger mean SST errors 
% than with a good offlux.  The offlux depends on the atmospheric forcing, 
% and is calculated from a previous run with that forcing.
%
% I'm using offlux for from Liz Douglass run ARCcb0.08 
%
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close

R = 'ARCc0.08';
E = '122';


%fnm_in  = 'offlux_237i+239i_ARCb0.08_11';
%fnm_in  = 'offlux_237i+239i_ARCb0.08_11';
fnm_in  = 'offlux_070+071i+100+102_2002-2012_mean_11';
fnm_out = 'offlux_070+071i+100+102_2002-2012_ZERO_11';

%ntopo1=09;
ntopo1=11; % topo for ARCc0.08
ntopo2=11; % topo for ARCc0.04
TV = sprintf('%2.2i',ntopo2);

pthin     = '/Net/kronos/ddmitry/hycom/ARCc0.08/offlux/';
pthout    = '/Net/kronos/ddmitry/hycom/ARCc0.08/offlux/';
%pthtopo   = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthtopo   = '/Net/kronos/ddmitry/hycom/ARCb0.08/topo/';
pthtopo04 = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/%s/topo_grid/',R);

% ARCb0.08 topo
flda=sprintf('%sdepth_ARCb0.08_11.a',pthtopo);
fldb=sprintf('%sdepth_ARCb0.08_11.b',pthtopo);
flga=sprintf('%sregional.grid.a',pthtopo);
flgb=sprintf('%sregional.grid.b',pthtopo);

% Read in topo from *a file:
f1=fopen(flgb);  % read I,J from regional.grid.b
aa=fgetl(f1);
dmm=aa(2:8);
IDM8=str2num(dmm);
aa=fgetl(f1);
dmm=aa(2:8);
JDM8=str2num(dmm);
IJDM8=IDM8*JDM8;
fclose(f1);

npad8=4096-mod(IJDM8,4096);

depth_fid=fopen(flda,'r');
%fseek(depth_fid,6*4*(npad+IJDM),-1) % this is confusing, should not be needed
[h,count]=fread(depth_fid,IJDM8,'float32','ieee-be');
%y=find(h>1e10);
%h(y)=nan;
HH8=reshape(h,IDM8,JDM8)';

fclose(depth_fid);

I=find(HH8>1e10);
HH8=-HH8;
HH8(I)=100;


 % read lon/lat from regional grid file
grid_fid=fopen(flga,'r');

[plon,count]=fread(grid_fid,IJDM8,'float32','ieee-be');
stat=fseek(grid_fid,4*(npad8+IJDM8),-1);
if stat<0
  error('Reading grid file ...');
end
[plat,count]=fread(grid_fid,IJDM8,'float32','ieee-be');

disp('Reading lat/lon  ...')
LON8=(reshape(plon,IDM8,JDM8))';
LAT8=(reshape(plat,IDM8,JDM8))';

fclose(grid_fid);

% --------------------------
I=find(LON8>180);
LON8(I)=LON8(I)-360;
I=find(LON8<-180);
LON8(I)=LON8(I)+360;


%fltopo = sprintf('%sdepth_ARCc0.08_%2.2i.nc',pthtopo,ntopo1);
%HHo   = nc_varget(fltopo,'Bathymetry');
%LATo  = nc_varget(fltopo,'Latitude');
%LONo = nc_varget(fltopo,'Longitude');
%[JDMo,IDMo]=size(HHo);
%IJDMo=IDMo*JDMo;
%npado=4096-mod(IJDMo,4096);

% New topo
fltopo_new=sprintf('%sdepth_%s_%s.nc',pthtopo04,R,TV);
HH   = nc_varget(fltopo_new,'Bathymetry');
LAT = nc_varget(fltopo_new,'Latitude');
LON = nc_varget(fltopo_new,'Longitude');
[mm,nn]= size(HH);

JDM = mm;
IDM = nn;
IJDM=IDM*JDM;
npad=4096-mod(IJDM,4096);
toto=ones(npad,1);

fprintf('%s domain, Topo=%s, ID=%i JD=%i\n',R,TV,IDM,JDM);

II8=[1:IDM]; % 
JJ8=[1:JDM];
II=[1:IDM];  
JJ=[1:JDM];

    
% ---------------------------------
%   Read offlux data:
% ---------------------------------
fina = sprintf('%s%s.a',pthin,fnm_in);
finb = sprintf('%s%s.b',pthin,fnm_in);
fida = fopen(fina,'r','ieee-be');
%fidb = fopen(finb,'r');
fprintf('Input file: %s\n %s\n',fina,finb);


frewind(fida);
A8 = fread(fida,IJDM8,'float32','ieee-be');
A8 = (reshape(A8,IDM8,JDM8))';
fclose(fida);

%keyboard

lm8=HH8*0;
lm8(HH8<0)=1;
lma=A8*0;
lma(A8<1e10)=1;
% Fill land with 0 for interpolation
Ia=find(A8>1e10);
A8(Ia)=0;

f_plt=0;
if f_plt==1
  cl1 = colormap_red(200);
  cl2 = flipud(colormap_cold(200));
  cmp = [cl2;cl1];

  A8(Ia)=nan;
  figure(1); clf;
  pcolor(A8); shading flat;
  colormap(cmp);
  caxis([-50 50]);
  colorbar

  ii1 = max(strfind(fina,'/'));
  title(fina(ii1+1:end),'Interpreter','none');
  btx = 'interp_offlux08to04.m';
  bottom_text(btx,'pwd',1); 
  
end

%keyboard

A = A8*0;

a1=min(min(A8));
a2=max(max(A8));
b1=min(min(A));
b2=max(max(A));
fprintf('  ARCc0.08 min/max: %8.5f %8.5f\n',a1,a2);
fprintf('  ARCc0.04 min/max: %8.5f %8.5f\n',b1,b2);

hg=2^100;
I4=find(HH>=0);
A(I4)=hg;
F=A';
F=reshape(F,IJDM,1);

fouta = sprintf('%s%s.a',pthout,fnm_out);
fouda = fopen(fouta,'w','ieee-be');
fwrite(fouda,F,'float32','ieee-be');
fwrite(fouda,toto,'float32','ieee-be');
fclose(fouda);
fprintf('Written files: %s\n',fouta);


foutb = sprintf('%s%s.b',pthout,fnm_out);
foudb = fopen(foutb,'wt');

aa = 'Heat Flux Offset (W/m^2)';
fprintf(foudb,[aa,'\n']);
aa = '-45.0*(expt_23.7+23.9_SSTann-PF_SSTann)';
fprintf(foudb,[aa,'\n']);
aa = '+150.0*(expt_23.7+23.9_ICEann-NOAA_ICEann)+150.0*(239_ICEann-NOAA_ICEann)';
fprintf(foudb,[aa,'\n']);
aa = 'anomaly between -2.0 and 2.0 degC';
fprintf(foudb,[aa,'\n']);
aa = sprintf('i/jdm =  %i  %i',IDM,JDM);
fprintf(foudb,[aa,'\n']);
aa = sprintf('min, max =  %10.7E    %10.7E',b1,b2);
fprintf(foudb,[aa,'\n']);
fclose(foudb);

fprintf('Written files: %s\n',foutb);







