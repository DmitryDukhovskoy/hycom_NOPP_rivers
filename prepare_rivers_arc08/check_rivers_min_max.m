%Check agreement of the min/max river runoff in *a, *b files
% Units:  m/s
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
startup;

close all
clear

PTH.data1='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/force/rivers/';
PTH.data = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/force/riversNCAR/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
%fltopo = sprintf('%sdepth_ARCc0.04_17DD.nc',pthtopo);
fltopo=sprintf('%sdepth_ARCc0.08_11.nc',pthtopo);


% Read annual fields at a time
yr1=2011;
yr2=2016;


% Get topo and grid:
HH   = nc_varget(fltopo,'Bathymetry');
LAT = nc_varget(fltopo,'Latitude');
LON = nc_varget(fltopo,'Longitude');
[m,n]= size(HH);
[mm,nn]= size(HH);
[DX,DY]=sub_dx_dy(LON,LAT);
ACell=DX.*DY;


% Reading rivers & bathymetry:
IDM=n;
JDM=m;
IJDM=IDM*JDM;
npad=4096-mod(IJDM,4096);

for YY=yr1:yr2;
%flriva=sprintf('%srivers_09.a',PTH.data);
%flrivb=sprintf('%srivers_09.b',PTH.data);
%flriva=sprintf('%srivers_11.a',PTH.data1);  % expt 11.0 clim. rivers, no Gr.
%flrivb=sprintf('%srivers_11.b',PTH.data1);
  flriva=sprintf('%srivers_11_NCAR_Gr_%4.4i.a',PTH.data,YY); % NCAR riv+Green.
  flrivb=sprintf('%srivers_11_NCAR_Gr_%4.4i.b',PTH.data,YY); % NCAR riv+Green.
  %fltopo=sprintf('%sdepth_ARCc0.08_09.nc',PTH.topo);

  fprintf('Reading %s\n',flriva);
  % Rivers:
  fida=fopen(flriva,'r');
  fidb=fopen(flrivb,'rt');

  for ii=1:4
    aa=fgetl(fidb);
    fprintf('%s\n',aa);
  end
  aa=fgetl(fidb);
  Is=strfind(aa,'=');
  S = sscanf(aa(Is+1:end),'%g');
  idm=S(1);
  jdm=S(2);
  if idm~=IDM | jdm~=JDM
    fprintf('IDM and JDM do not agree with topo: %i, %i\n',idm,jdm);
  end
  
  for MM=1:12  
    fprintf('Year %i, Mo=%i\n',YY,MM);

    aa=fgetl(fidb);
    if isempty(aa) | aa==-1; fprintf('--- EoF -----\n'); break; end;
    Is=strfind(aa,'=');
    fld=aa(1:Is-2);
    S = sscanf(aa(Is+1:end),'%g');
    imo  = S(1);
    smin = S(2);
    smax = S(3);

    fseek(fida,(MM-1)*4*(npad+IJDM),-1);
    [A,counta]=fread(fida,IJDM,'float32','ieee-be');
    toto = fread(fida,npad,'float32','ieee-be');
    I=find(A<1e20); % although in river files no huge values
    amin=min(A(I));
    amax=max(A(I));

    fprintf('%i/%i .a, .b min = %16.7e %16.7e\n',...
	  YY,MM,amin,smin);
    fprintf('%i/%i .a, .b max = %16.7e %16.7e\n',...
	  YY,MM,amax,smax);

    dmn=abs(1-amin/smin);
    dmx=abs(1-amax/smax);
%keyboard
    if dmn>1e-4
      fprintf(' *****   ERR:  .a and .b not consistent dmn=%9.7f\n',dmn);
      keyboard;
    end
    if dmx>0.0001
      fprintf(' *****  ERR: .a and .b not consistent dmx=%9.7f\n',dmx);
      keyboard;
    end

  end; % months
  fclose(fida);
  fclose(fidb);
end

fprintf('Min/Max match, all done ...\n');







