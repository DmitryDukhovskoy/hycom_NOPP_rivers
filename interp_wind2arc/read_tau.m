% Read/plot  wind stress forcing for checking
% Wind stress fields NRL CFSR
% 1-hr fields
% sea mask over land
% QuikSCAT corrected wind speeds
% wind stress is computed by NRL
% using bulk formulae
%
% tau*_yyyA.a - Years/months are in hycom notation !!!
% Time is NRL HYCOM time: days since dec 31 1900
% Direct access Fortran file
%
% Dmitry Dukhovskoy COAPS FSU March 2016
%     October 2019: few changes for Theresa
%                   
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

MNTH={'a','b','c','d','e','f','g','h','i','j','k','l'};

% Specify time to read/plot
dplot = datenum(2005,01,31,23,0,0);
DV=datevec(dplot);

fprintf('Reading wind stress %s\n',datestr(dplot));

YR=DV(1);
yy=sprintf('%3.3i',YR-1900);
mnth=MNTH{DV(2)};

[dS,dE,dC] = get_dtime(dplot);

tread=dC; % HYCOM time

pthin    = sprintf('/nexsan/archive/ARCc0.08_110/data_natl/%4.4i/',YR);
PTH.topo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
PTH.data = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/data_ice_restart/';

fltopo=sprintf('%sdepth_ARCc0.08_11.nc',PTH.topo);
HH   = nc_varget(fltopo,'Bathymetry');
alat = nc_varget(fltopo,'Latitude');
elon = nc_varget(fltopo,'Longitude');
LAT  = alat;
LON  = elon;
[m,n]= size(HH);
[mm,nn]= size(HH);


%ID=1600;
%JD=2520;

% Find record # to read:
fina =sprintf('%stauewd_%s%s.a',pthin,yy,mnth);
fina2=sprintf('%staunwd_%s%s.a',pthin,yy,mnth);
finb =sprintf('%stauewd_%s%s.b',pthin,yy,mnth);

f1 = fopen(finb,'r');  % read I,J from *.b
for nl=1:4
  aa=fgetl(f1);
%  disp(aa);
end

aa=fgetl(f1);
[c1,c2,ID, JD] = strread(aa,'%s%s%d%d');
IJDM=ID*JD;
npad=4096-mod(IJDM,4096);
toto=ones(npad,1);
fprintf('Retrieved dimensions: I=%i, J=%i\n',ID,JD);

irec=0;
irec0=0;
while ~irec0
  aa=fgetl(f1);
  if ~char(aa); break; end; % E-O-F
  ir=strfind(aa,'=');
  bb=aa(ir+1:end);
  irec=irec+1;
  dmm = textscan(bb,'%f%f%f%f');
  dtime=dmm{1};
  if irec==1, 
    dtime1=dtime; 
  end;
  if abs(dtime-tread)<1e-4,
    fprintf('Found Rec. # %i\n',irec);
    irec0=irec;
    break;
  end
end

if irec0==0
  fprintf('Could not find time tread=%9.5f\n',tread);
  fprintf('Time range in file: %9.5f - %9.5f\n',dtime1,dtime);
  error('**** ERROR:  Check tread  *****');
end


fclose(f1);

%keyboard

fid1=fopen(fina,'r');
frewind(fid1);
k0=irec0-1;
stat=fseek(fid1,k0*(IJDM+npad)*4,-1);
dmm=fread(fid1,IJDM,'float32','ieee-be');  % read 2D field (1 layer)
dmm=reshape(dmm,ID,JD);
fclose(fid1);

Fx=dmm';

fid1=fopen(fina2,'r');
frewind(fid1);
k0=irec0-1;
stat=fseek(fid1,k0*(IJDM+npad)*4,-1);
dmm=fread(fid1,IJDM,'float32','ieee-be');  % read 2D field (1 layer)
dmm=reshape(dmm,ID,JD);
fclose(fid1);

Fy=dmm';



figure(1); clf;
pcolor(Fx); shading flat;
hold on;
contour(HH,[0 0],'k');

colorbar
%tts=sprintf('Day %9.5f, %s',dtime,fina);
tts=sprintf('taux, Day %9.5f, %s',dtime,datestr(dplot));
title(tts,'Interpreter','none');

figure(2); clf;
pcolor(Fy); shading flat;
hold on;
contour(HH,[0 0],'k');

colorbar
%tts=sprintf('Day %9.5f, %s',dtime,fina2);
tts=sprintf('tauy, Day %9.5f, %s',dtime,datestr(dplot));
title(tts,'Interpreter','none');






