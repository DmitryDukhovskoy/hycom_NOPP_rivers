% Read hycom restart fields - initial fields
% for ARCc0.04 GOFS3.5
%
% One set of initial fields is prepared from my old AO HYCOM-CICE0.04 experiment 01.2
% another is from the Global 0.04 GLBc HYCOM-CICE5, GOFS3.5 
%
% Read all fields and check min/max with *b
% Check Land mask for sal temp
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
%addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

hg = 1e20;

%frst = '0.04 AO HYCOM-CICE4-012';
frst = 'GLBc0.04-GOFS3.5';


%pthbin  = '/Net/mars/ddmitry/hycom/ARCc0.04/restart/';
%pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/topo_grid/';
%fltopo = sprintf('%sdepth_ARCc0.04_17DD.nc',pthtopo);
if strncmp(frst,'GLBc0.04',8)
  pthbin  = '/Net/kronos/ddmitry/hycom/ARCc0.04/restart_022/';  % expt 02.2 - HYCOM-CICE5, GOFS3.5
else
  pthbin  = '/Net/kronos/ddmitry/hycom/ARCc0.04/restart_012/';  % expt 01.2 - HYCOM-CICE4, GOFS3.1
end
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/topo_grid/';
fltopo = sprintf('%sdepth_ARCc0.04_17DD.nc',pthtopo);

fprintf('Reading topo: %s\n',fltopo);
HH  = nc_varget(fltopo,'Bathymetry');
LAT = nc_varget(fltopo,'Latitude');
LON = nc_varget(fltopo,'Longitude');
[mm,nn]= size(HH);
[m,n]= size(HH);
IDM=n;
JDM=m;
IJDM=IDM*JDM;
npad=4096-mod(IJDM,4096);
toto=ones(npad,1);
Ib = find(HH<0);

if strncmp(frst,'GLBc0.04',8)
  fina = sprintf('%sGLBc2ARCc004_restart_117a.a',pthbin);
  finb = sprintf('%sGLBc2ARCc004_restart_117a.b',pthbin);
else
  fina = sprintf('%srestart_from012_008a_116a.a',pthbin);
  finb = sprintf('%srestart_from012_008a_116a.b',pthbin);
end

fida = fopen(fina,'r');
fidb = fopen(finb,'rt');

aa=fgetl(fidb);
fprintf('%s\n',aa);
aa=fgetl(fidb);
fprintf('%s\n',aa);

% Start reading fields:
while ~isempty(aa),
  aa=fgetl(fidb);
  if isempty(aa) | aa==-1; fprintf('--- EoF -----\n'); break; end;

  Is=strfind(aa,'=');
  fld=aa(1:Is-2);
  S = sscanf(aa(Is+1:end),'%g');
  smin = S(3);
  smax = S(4);
  Lr   = S(1);
  TLev = S(2);
  
  [A,counta]=fread(fida,IJDM,'float32','ieee-be');
  toto = fread(fida,npad,'float32','ieee-be');
  I=find(A<0.1*hg);
  amin=min(A(I));
  amax=max(A(I));
  
  fprintf('\n %s  Layer %i, Time %i, .a, .b min = %16.7e %16.7e\n',...
	  fld,Lr,TLev,amin,smin);
  fprintf(' %s  Layer %i, Time %i, .a, .b max = %16.7e %16.7e\n',...
	  fld,Lr,TLev,amax,smax);
  
  dmn=abs(1-amin/smin);
  dmx=abs(1-amax/smax);
  
  if dmn>1e-4
    fprintf('.a and .b not consistent dmn=%9.7f\n',dmn);
    keyboard;
  end
  if dmx>0.0001
    fprintf('.a and .b not consistent dmx=%9.7f\n',dmx);
    keyboard;
  end
  
% Check nan values in ocean points:
  B = reshape(A,size(HH'))';
  Inan = find(isnan(B(Ib)) | B(Ib)>=0.1*hg);
  if ~isempty(Inan),
    fprintf('Found %i nan/huge values in ocean points...\n',length(Inan));
    keyboard;
  end
  

end;

fclose(fida);
fclose(fidb);

fprintf('No problems found in *a and *b min/max values\n')










