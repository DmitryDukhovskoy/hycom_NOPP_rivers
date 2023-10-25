% tbaric - thermobaric reference state
% To generate tbaric fields
% use tbaric.com
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close

R = 'ARCc0.04';
%R = 'ARCc0.08';
% E = '110'; % arc08 control run 1993-2014 with tracers
E = '010'; % ARCc0.04 expt - test run 01.0

ntopo1=09;  % topo used in old ARCc0.08 experiments with 32 layers
ntopo2=11;  % topo used in ARCc0.08 experiments with 41 layers
ntopo3=17;  % topo used in ARCc0.04

switch (R),
 case('ARCc0.08');
  TV = sprintf('%2.2i',ntopo2);
 case('ARCc0.04');
  TV = sprintf('%2.2iDD',ntopo3); % corrected - open strait in CAA
end


ptharc  = sprintf('/Net/kronos/ddmitry/hycom/%s/tmp_files/',R);
%pthglb  = '/Net/kronos/ddmitry/hycom/GLBb0.08/expt_69.1/restart/';
pthin  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/%s/relax/010/',R); % relax 32lyrs
pthout = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/%s/relax/%s/41layers_T%s/',R,E,TV);
pthtopo = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/%s/topo_grid/',R);


% Get new topo and grid:
fltopo_new=sprintf('%sdepth_%s_%s.nc',pthtopo,R,TV);
HH   = nc_varget(fltopo_new,'Bathymetry');
alat = nc_varget(fltopo_new,'Latitude');
elon = nc_varget(fltopo_new,'Longitude');
LAT  = alat;
LON  = elon;
[mm,nn]= size(HH);
% Old topo:
%fltopo_old=sprintf('%sdepth_ARCc0.08_%2.2i.nc',pthtopo,ntopo1);
%HHo  = nc_varget(fltopo_old,'Bathymetry');



% Read in sigma layers:
%fina = sprintf('%stbaric.a',pthin);
%finb = sprintf('%stbaric.b',pthin);
fina = sprintf('%stbaric_T%s.a',pthout,TV);
finb = sprintf('%stbaric_T%s.b',pthout,TV);

%fida = fopen(fina,'r','ieee-be');
fida = fopen(fina,'r','ieee-be');
fidb = fopen(finb,'r');
aa = fgetl(fidb);
aa = fgetl(fidb);
ie = strfind(aa,'=')+1;
CC = textscan(aa(ie:end),'%u%u');
fclose(fidb);

IDM = CC{1};
JDM = CC{2};
IJDM=IDM*JDM;
npad=4096-mod(IJDM,4096);
toto=ones(npad,1);

fprintf('ARCc domain, ID=%i JD=%i\n',IDM,JDM);

% Read *a file:
fida = fopen(fina,'r');
A = fread(fida,IJDM,'float32','ieee-be');
A = (reshape(A,IDM,JDM))';
fclose(fida);

A(A>1e6)=nan;

figure(2); clf;
%pcolor(A); shading flat;
%colorbar

contour(HH,[0 0],'k');
hold on;
contour(A,[1 1],'r');
contour(A,[2 2],'g');


