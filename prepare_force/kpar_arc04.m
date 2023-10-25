% Interpolate SeaWifs ocean color
% kpar.a to ARCc0.04
% also check range in created fields by using Alan's code:
% /nexsan/people/ddmitry/Net_ocean/HYCOM/hycom/ALL4/bin/hycom_range kpar_arc04.a 3200 5040
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

pthin   = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/force/seawifs/';
pthout  = '/Net/mars/ddmitry/hycom/ARCc0.04/force/';
pthtopo2= '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';
fltopo2 = sprintf('%sdepth_ARCc0.04_17DD.nc',pthtopo2); % 

IDM1=1600;
JDM1=2520;
IJDM1=IDM1*JDM1;
npad1=4096-mod(IJDM1,4096);
toto1=ones(npad1,1);

% Get new topo and grid:
HH2  = nc_varget(fltopo2,'Bathymetry');
LAT2 = nc_varget(fltopo2,'Latitude');
LON2 = nc_varget(fltopo2,'Longitude');
[m2,n2]= size(HH2);
IDM2=n2;
JDM2=m2;
IJDM2=IDM2*JDM2;
npad2=4096-mod(IJDM2,4096);
toto2=ones(npad2,1);

% For interpolation:
IIo=[1:2:IDM2+2]; % extend outside the new domain
JJo=[1:2:JDM2+2];
II=[1:IDM2];
JJ=[1:JDM2];


fina  = sprintf('%skpar_fixed.a',pthin);
finb  = sprintf('%skpar.b',pthin);
fouta = sprintf('%skpar_arc04.a',pthout);
foutb = sprintf('%skpar_arc04.b',pthout);
f1    = fopen(finb,'r');  % read I,J from *.b
fidaw = fopen(fouta,'w');
fidbw = fopen(foutb,'wt');

fprintf('Reading fields: %s\n',fina);
fprintf('Writing fields to: %s\n',fouta);

for nl=1:4
  aa=fgetl(f1);
  disp(aa);
  fprintf(fidbw,[aa,'\n']);
end

aa=fgetl(f1);
%[c1,c2,ID1, JD1] = strread(aa,'%s%s%d%d');
aa=sprintf('i/jdm = %i %i',IDM2,JDM2);
fprintf(fidbw,[aa,'\n']);

fclose(f1);

fid1=fopen(fina,'r');
frewind(fid1);
for ii=1:12
  fprintf('Interpolating, ii = %i\n',ii);
  dmm=fread(fid1,IJDM1,'float32','ieee-be');  % read 2D field (1 layer)
  dm1=fread(fid1,npad1,'float32','ieee-be');  % read npad 
  amin=min(min(dmm));
  amax=max(max(dmm));
  fprintf('  Old kpar: ii=%i, min=%13.10e, max=%13.10e\n',ii,amin,amax);
  dmm=reshape(dmm,IDM1,JDM1);
  F1=dmm';
  
  if isempty(dmm) | isempty(dm1),
    error('EOF - input file');
  end
  
  F1(JDM1+1,:) = F1(JDM1,:);
  F1(:,IDM1+1) = F1(:,IDM1);
  F2 = interp2(IIo,JJo',F1,II,JJ');
  Inan = find(isnan(F2));
  if ~isempty(Inan),
    error('There are nan"s in the interpolated kpar');
  end
  
  A  = F2';
  A  = reshape(A,IJDM2,1);
  fwrite(fidaw,A,'float32','ieee-be');
  fwrite(fidaw,toto2,'float32','ieee-be');
  amin=min(min(A));
  amax=max(max(A));
  fprintf('  New kpar: ii=%i, min=%13.10e, max=%13.10e\n',ii,amin,amax);
  aa = sprintf('    kpar: month,range = %2.2i  %10.7E  %10.7E',ii,amin,amax);
  fprintf(fidbw,[aa,'\n']);
  
%  keyboard
end;

fprintf(' Done \n');
fclose(fidaw);
fclose(fidbw);
fclose(fid1);
%exit