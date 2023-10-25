% Analyze GIN Seas climatology
% from NOAA 
% based on observations
% Greenland-Iceland-Norwegian Seas Regional Climatology version 2 (GINS RC v2)
% https://www.nodc.noaa.gov/OC5/regional_climate/gin-seas-climate/about_gin.html
% The GINS RC v2 is a set of six decadal mean fields for temperature and 
% salinity within the 50°N and 85°N and 45°W to 15°E domain comprising 
% the Greenland, Iceland, and Norwegian Seas and adjacent 
% areas of the northern North Atlantic and Arctic Ocean. 
% The six decadal climatologies of temperature and salinity use 
% all data available in WOD13 from 1955 to 2012. 
% There are five 10-year 
% (1955-1964, 1965-1974, 1975-1984, 1985-1994, and 1995-2004) 
% climatologies and an 8-year (2005-2012) climatology 
% that were calculated for annual, seasonal, 
% and monthly time periods. 
% Seasons are as follows: Winter (January-March), Spring (April-June), 
% Summer (July-September), and Fall (October-December). 
% Additionally, there is a 58-year climatology defined 
% as "all averaged decades," which was computed by 
% averaging all six decadal climatologies and 
% represents the long-term mean state of the ocean in the GINS region.
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers
startup;

% annual statistical mean S, by decades
% 2005-2012
url_file = ['https://data.nodc.noaa.gov/thredds/',...
	    'dodsC/nodc/archive/data/0170742/DATA/',...
	    'salinity/netcdf/A5B2/0.10/gins_A5B2_s00_10.nc'];
%ncdisp(url_file);

X = ncread(url_file,'lon');
Y = ncread(url_file,'lat');
ZZ= ncread(url_file,'depth');
ZZ= -abs(ZZ);
nz= length(ZZ);
%S = ncread(url_file,'s_mn'); % mean of unflagged interpolated values
S = ncread(url_file,'s_an'); % objectively analyzed s fields
S = shiftdim(S,2);

A=[];
for k=1:nz
  dmm=squeeze(S(k,:,:));
  A(k,:,:)=dmm';
end
S = A;

LR=[0,-50;...
    -50,-150;...
    -150,-300];

dZ = abs(diff(ZZ));

for ilv=1:3
  z1=LR(ilv,1);
  z2=LR(ilv,2);
  k1=min(find(ZZ<=z1));
  k2=max(find(ZZ>=z2));
  
  dZm=0;
  [ll,mm,nn] = size(S);
  Sm=zeros(mm,nn);
  for kk=k1:k2
    ss = squeeze(S(kk,:,:));
    Sm=Sm+ss*dZ(kk);
    dZm=dZm+dZ(kk);
  end
  Sm=Sm./dZm;
end


