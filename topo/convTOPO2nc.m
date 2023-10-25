% for visualization purposes, convet
% HYCOM bathymetry to netCDF
% is subregion is extracted form HYCOM output, then
% need to save TOPO subregion as well
%

RG='ARCc0.08';
%RG='ARCb0.08';
%T=11;            % topo #
T=7;            % topo #
EX='110';       % experiment #
pthT = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/%s/topo_grid/',RG);
%pthout=sprintf('/Net/bimonthly/ddmitry/%s/%s/topo_nc/',RG,EX);
pthout=pthT;

% Topo/grid files, input:
flda=sprintf('%sdepth_%s_%2.2i.a',pthT,RG,T);
flga=sprintf('%sregional.grid.a',pthT);
flgb=sprintf('%sregional.grid.b',pthT);

% Output file:
%fnmB = sprintf('depth_%s_%2.2i_Nord',RG,T);
fnmB = sprintf('depth_%s_%2.2i',RG,T);
fnm=[pthout,fnmB];
fcdl = [fnm,'.cdl'];


% -------------------------------------------------
%       SUBREGIONING:
% -------------------------------------------------
%js=401;  % start index j
%je=1700;  % end index j
%is=601;
%ie=1600;
js=1;
je=2520;
is=1;
ie=1600;

nn=ie-is+1;
mm=je-js+1;

%if ~exist(pthout,'dir');
%  str=sprintf('! mkdir -v %s ',ftout);
%  eval(str);
%end;

fprintf('Region: %s, Experiment: %s, Topo: %2.2i\n',RG,EX,T);
fprintf('Input Topo file: %s \n',flda);
fprintf('Start i=%i, end i=%i, start j=%i, end j=%i,  n=%i, m=%i\n',is,ie,js,je,nn,mm);


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
%keyboard

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
y=find(h>1e10);
h(y)=nan;
h=reshape(h,IDM,JDM)';

fclose(depth_fid);

%keyboard
% ==========================
%    subregion domain:
% ==========================
LON=plon(js:je,is:ie);
LAT=plat(js:je,is:ie);
HH =h(js:je,is:ie);
if min(min(HH))>0;
  HH=-HH;
end;
I=find(isnan(HH));
HH(I)=100;

%========================
% Create cdl file:

fid = fopen(fcdl,'wt');
fprintf(fid,'%s \n',sprintf('netcdf %s {',fnmB));
fprintf(fid,'\n');
fprintf(fid,'%s \n','dimensions:');
fprintf(fid,'%s \n',sprintf('    Y = %i, X = %i ;',mm,nn));
fprintf(fid,'\n');
fprintf(fid,'%s \n','variables:');
fprintf(fid,'%s \n','    float Bathymetry(Y,X), Latitude(Y,X), Longitude(Y,X) ;');
fprintf(fid,'\n');
% Variable attributes:
fprintf(fid,'%s \n','    Bathymetry:standard_name = "sea floor topo" ;');
fprintf(fid,'%s \n','    Bathymetry:units = "m" ;');
fprintf(fid,'%s \n','    Bathymetry:positive = "above sea surface" ;');
fprintf(fid,'%s \n','    Bathymetry:land = "100" ;');
fprintf(fid,'%s \n','    Latitude:standard_name = "latitude" ;');
fprintf(fid,'%s \n','    Latitude:units = "degrees_north" ;');
fprintf(fid,'%s \n','    Longitude:standard_name = "longitude" ;');
fprintf(fid,'%s \n','    Longitude:units = "degrees_east" ;');
fprintf(fid,'\n\n');
fprintf(fid,'%s \n','//global attributes:');
fprintf(fid,'%s \n','                :Conventions = "CF-1.0" ;');
fprintf(fid,'%s \n','                :title = "HYCOM ARCc0.08 FSU COAPS" ;');
fprintf(fid,'%s \n','                :reference1 = "NRL SSC A. Wallcraft, P. Posey" ;');
fprintf(fid,'%s \n','                :reference2 = "COAPS FSU D.Dukhovskoy" ;');
fprintf(fid,'%s \n',sprintf('                :topography = "T%2.2i" ;',T));
fprintf(fid,'%s \n',sprintf('                :experiment = "%s" ;',EX));
fprintf(fid,'%s \n',sprintf('                :created = "%s" ;',date));
fprintf(fid,'}');
fclose(fid);

stt=sprintf('!ncgen -b %s',fcdl);
eval(stt);
st2=sprintf('! mv %s.nc %s.',fnmB,pthout);
eval(st2);

fnc = [pthout,fnmB,'.nc'];
disp('Writing netCDF topo file ...');
disp(fnc);
nc_varput(fnc,'Bathymetry',HH);
nc_varput(fnc,'Latitude',LAT);
nc_varput(fnc,'Longitude',LON);

disp('DONE ');






