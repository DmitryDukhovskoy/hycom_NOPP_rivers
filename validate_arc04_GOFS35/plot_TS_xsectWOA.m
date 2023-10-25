% Vertical section of T/S from World Ocean Atlas 2018 1/4 degree
% Monthly
% File naming conventions: [V][TT][FF][GG].[EXT]
%where:
%[V] - variable
%[TT] - time period
%[FF] - field type
%[GG] - grid (01- 1°, 04 - 1/4° 10 - 1/10°)
%[EXT] - file extention
%Note: '.dat' - ASCII; '.csv' - comma separated value; '.dbf',
% '.shp', '.shx' - ArcGIS shape files; '.nc' - netCDF files

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

%
imo = 7;
pfld = 'temp';
fld0 = pfld;

pT='https://www.ncei.noaa.gov/thredds-ocean/dodsC/ncei/woa/temperature/decav/0.25/';
pS='https://www.ncei.noaa.gov/thredds-ocean/dodsC/ncei/woa/salinity/A5B7/0.25/';

if strncmp(pfld,'temp',4)
  flnm=sprintf('woa18_decav_t%2.2i_04.nc',imo);
  pth=pT;
else
%woa18_A5B7_s07_04.nc
  flnm=sprintf('woa18_decav_s%2.2i_04.nc',imo);
  pth=pS;
end

finp=sprintf('%s%s',pth,flnm);



