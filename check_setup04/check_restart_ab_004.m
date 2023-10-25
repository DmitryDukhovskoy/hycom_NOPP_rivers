% Read hycom restart fields
% for ARCc0.04
% Read all fields and check min/max with *b
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

hg = 1e20;

%pthbin  = '/Net/mars/ddmitry/hycom/ARCc0.04/restart/';
%pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';
%fltopo = sprintf('%sdepth_ARCc0.04_17DD.nc',pthtopo);
pthbin  = '/Net/mars/ddmitry/hycom/ARCc0.04/restart/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';
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

fina = sprintf('%srestart_105a.a',pthbin);
finb = sprintf('%srestart_105a.b',pthbin);

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










