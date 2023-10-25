% Read sea ice monthly area for the Arctic Ocean
% bootstrap form NSIDC Earthdata website
% Daily Arctic and Antarctic ice-covered area and total sea ice extent
%  ASCII text data files list the total ice-covered area (km2) and total
%  sea ice extent (km2) for both the Bootstrap Sea Ice Concentrations 
% from Nimbus-7 SMMR and DMSP SSM/I-SSMIS and/or the NASA Team Sea Ice 
% Concentrations from Nimbus-7 SMMR and DMSP SSM/I-SSMIS Passive Microwave Data algorithms.
% TSIE and TSIA from NSIDC
% https://nsidc.org/data/NSIDC-0192/versions/3/print
%
% data at 
%
% https://daacdata.apps.nsidc.org/pub/DATASETS/nsidc0192_seaice_trends_climo_v3
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

pthiconc   = '/nexsan/people/ddmitry/Net_ocean/SeaIce/';
fice = sprintf('%sgsfc.bootstrap.month.area.1978-2020.n.txt',pthiconc);
fprintf('Reading %s\n',fice);

fid = fopen(fice,'r');
frewind(fid);
C = textscan(fid,'%s','HeaderLines',1);
fclose(fid);

% Columns:
% Year  Mon  Ver  TotalArc  Okhotsk   Bering   Hudson   Baffin  Grnland  BarKara  ArctOcn  CanArch   StLawr
Nflds = 13;
Nrec = length(C{1})/Nflds;
icnt=0;
for ii=1:Nrec
  iyr  = (ii-1)*Nflds+1;
  imo  = (ii-1)*Nflds+2;
  itot = (ii-1)*Nflds+4;
  iarc = (ii-1)*Nflds+11;

  Year    = str2double(C{1}(iyr));
  mo      = str2double(C{1}(imo));
  TotArc  = str2double(C{1}(itot));
	ArctOcn = str2double(C{1}(iarc)); 

  icnt = icnt+1;
  ICE.Year(icnt,1)            = Year;
  ICE.Month(icnt,1)           = mo;
  ICE.TotArctArea_km2(icnt,1) = TotArc;
  ICE.ArctOcnArea_km2(icnt,1) = ArctOcn;

  fprintf('Year =%i Month=%2.2i TotArtArea =%8.4fx10^6 km2\n',Year,mo,TotArc*1e-6);
end

fmat = sprintf('%sBootstrap_Arctic_IceArea.mat',pthiconc);
fprintf('Saving %s\n',fmat);
save(fmat,'ICE');


