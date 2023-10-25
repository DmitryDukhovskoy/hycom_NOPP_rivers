% Prepare bath and grid infor
%
% Plot regino map
% with bath + Gr runoff
% + obs points
% Plot river sources
% Check onland river sources
% Units:  m/s
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
startup;

close all
clear

PTH.data1='/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/force/rivers/';
PTH.data = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/force/riversNCAR/';
PTH.topo='/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
PTH.mat='/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_theresa/';

expt = 112; % NCAR + Greenland runoff
%expt = 110; % river climatology
YY=2015;
MM = 7; % month
YYx=YY;
MMx=MM+1;
if MMx==13,
  MMx=1;
  YYx=YY+1;
end

mdays = datenum(YYx,MMx,1)-datenum(YY,MM,1);

%flriva=sprintf('%srivers_09.a',PTH.data);
%flrivb=sprintf('%srivers_09.b',PTH.data);
%flriva=sprintf('%srivers_11_Greenland_%i.a',PTH.data,YY);
%flrivb=sprintf('%srivers_11_Greenland_%i.b',PTH.data,YY);
switch(expt)
 case(110);
  flriva=sprintf('%srivers_11.a',PTH.data1);  % expt 11.0 clim. rivers, no Gr.
  flrivb=sprintf('%srivers_11.b',PTH.data1);
 case(112)
  flriva=sprintf('%srivers_11_NCAR_Gr_%4.4i.a',PTH.data,YY); % NCAR riv+Green.
  flrivb=sprintf('%srivers_11_NCAR_Gr_%4.4i.b',PTH.data,YY); % NCAR riv+Green.
end

%fltopo=sprintf('%sdepth_ARCc0.08_09.nc',PTH.topo);
fltopo=sprintf('%sdepth_ARCc0.08_11.nc',PTH.topo);

fprintf('Reading %s\n',flriva);
fprintf('Expt=%3.3i, Year %i, Mo=%i\n',expt,YY,MM);

% Get topo and grid:
HH   = nc_varget(fltopo,'Bathymetry');
LAT = nc_varget(fltopo,'Latitude');
LON = nc_varget(fltopo,'Longitude');
[m,n]= size(HH);
[mm,nn]= size(HH);
[DX,DY]=sub_dx_dy(LON,LAT);
ACell=DX.*DY;

% Get Greenland contour:
GC = sub_greenl_isobath(HH,LON,LAT);

GRID.Topo=HH;
GRID.Long=LON;
GRID.Latit=LAT;
GRID.DX_spacing=DX;
GRID.DY_spacing=DY;

fmat=sprintf('%shycom008_grid_GreenlContour',...
	     PTH.mat);
fprintf('Saving %s\n',fmat);
save(fmat,'GC','GRID');



