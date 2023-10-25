% Create cice atm. forcing fields
% Climatology precipitation
% Using ARCc0.08 cice.prec_lanl_12.r file
%  mm/mo

addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close

R = 'ARCc0.04';
E = '010';
%ntopo1=09;
ntopo1=11; % topo for ARCc0.08
ntopo2=17; % topo for ARCc0.04
TV = sprintf('%2.2iDD',ntopo2);

ptharc    = '/Net/kronos/ddmitry/hycom/ARCc0.08/tmp_files/';
%pthglb  = '/Net/kronos/ddmitry/hycom/GLBb0.08/expt_69.1/restart/';
pthin     = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/force/';
pthout    = sprintf('/Net/mars/ddmitry/hycom/%s/force/',R);
pthtopo   = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthtopo04 = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/%s/topo_grid/',R);

% Old topo:
fltopo = sprintf('%sdepth_ARCc0.08_%2.2i.nc',pthtopo,ntopo1);
HHo   = nc_varget(fltopo,'Bathymetry');
%LATo  = nc_varget(fltopo,'Latitude');
%LONo = nc_varget(fltopo,'Longitude');
[JDo,IDo]=size(HHo);

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

% arc0.08 grid is every second grid point of arc0.04
% arc0.04 last points are ~3.5 km outside the arc0.08 grid
IIo=[1:2:IDM+1]; % for interpolation, add extra col/row
JJo=[1:2:JDM+1];
II=[1:IDM];
JJ=[1:JDM];

% -------------------
%   Read precip data
% -------------------
fin = sprintf('%scice.prec_lanl_12.r',pthin);
fid = fopen(fin,'r','ieee-be');

fout = sprintf('%scice.prec_lanl_12_arc04.r',pthout);
foud = fopen(fout,'w','ieee-be');

frewind(fid);
for ii=1:12
  fprintf('Interpolating, ii = %i\n',ii);
%  dmm=fread(fid,1,'float64');
  dmm=fread(fid,[IDo,JDo],'float64');
  A08 = dmm';
% Add extra column/row for interpolation
  A08(:,end+1)=A08(:,end);
  A08(end+1,:)=A08(end,:);
%  dmm=fread(fid,1,'float64');
  if (isempty(dmm)); 
    fprintf('E-o-F, month = %i\n',ii);
    keyboard
  end
  A = interp2(IIo,JJo',A08,II,JJ');
  fwrite(foud,A','float64');

  inan=find(isnan(A));
  if ~isempty(inan)
    error('Found nans in interpolated rhoa');
  end

  a1=min(min(A08));
  a2=max(max(A08));
  b1=min(min(A));
  b2=max(max(A));
  fprintf('  ARCc0.08 min/max: %8.5f %8.5f\n',a1,a2);
  fprintf('  ARCc0.04 min/max: %8.5f %8.5f\n',b1,b2);
  

end

fclose(fid);
fclose(foud);

fprintf('Written files: %s\n',fout);







