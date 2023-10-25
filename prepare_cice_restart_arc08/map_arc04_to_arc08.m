% Remap CICE restart from GOFS3.5-type
% ARCc0.04 --> ARCc0.08

addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear
close all


pthnc = '/nexsan/people/ddmitry/hycom/ARCc0.04_022/incoming/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthtopo4 = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';

%narc04 = 'cice.restart04_117f.nc';
%narc08 = 'cice.restart08_117f_DD.nc';
narc04 = 'cice.restart_117d.nc';
narc08 = 'cice.restart08_117d_DD.nc';

IDM = 1600;
JDM = 2520;

%
% Topography:
ftopo = sprintf('%sdepth_ARCc0.08_11.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
Ilnd = find(HH>=0);

%
% 0.04 grid:
ftopo04 = sprintf('%sdepth_ARCc0.04_17DD.nc',pthtopo4); % 
HH4  = nc_varget(ftopo04,'Bathymetry');
LN4 = nc_varget(ftopo04,'Longitude');
LT4 = nc_varget(ftopo04,'Latitude');
[m4,n4] = size(LN4);

%
% Find ocean pnts in 0.08 that are lands in 0.04:
H48=HH4(1:2:m4,1:2:n4);
Ilo=find(HH<0 & H48>=0);  % land points in HYCOM0.04
[JJ,II]=ind2sub(size(HH),Ilo);

fin  = sprintf('%s%s',pthnc,narc04);
fout = sprintf('%s%s',pthnc,narc08);

info_input = ncinfo(fin);
info_input.Dimensions(1).Length = IDM;
info_input.Dimensions(2).Length = JDM;


%
% INPUT fields:
fprintf('Initial file: %s\n',fin);
finid = netcdf.open(fin,'NC_NOWRITE');
[ndims,nvars,ngatts,unlimdimid] = netcdf.inq(finid);
for iv=1:nvars
  [varname vartype vardimIDs varatts] = netcdf.inqVar(finid,iv-1);
  fprintf('%i %s\n',iv,varname);
  BI(iv).Name = varname;
  BI(iv).Ndim = length(vardimIDs);
end
netcdf.close(finid);

%
% Changes dimensions:
for ivar = 1:nvars
  nd = length(info_input.Variables(ivar).Dimensions);
  
  ini = [];
  inj = [];
  for idim=1:nd
    dnm = info_input.Variables(ivar).Dimensions(idim).Name;
    if strncmp(dnm,'ni',2)
      ini=idim;
      info_input.Variables(ivar).Dimensions(idim).Length=IDM;
    elseif strncmp(dnm,'nj',2)
      inj=idim;
      info_input.Variables(ivar).Dimensions(idim).Length=JDM;
    end
  end
  info_input.Variables(ivar).Size(ini) = IDM;
  info_input.Variables(ivar).Size(inj) = JDM;
end


% Writing output netcdf
fprintf('Creating %s\n',fout);
if exist(fout,'file')
  scmd=sprintf('/bin/rm %s',fout);
  system(scmd);
end

ncwriteschema(fout,info_input);

for ivar = 1:nvars
  nname = BI(ivar).Name;

% Find the closest point:
% which is every 2nd points on ARCc0.04
  fprintf('Writing %s\n',nname);
  Ain = nc_varget(fin,nname);
  ndm = length(size(Ain));
  A = [];
  if ndm==2
    A = Ain(1:2:m4,1:2:n4);
    A(Ilo)  = 0.0;
    A(Ilnd) = 0.0;
  else
    [ncc,n1,n2] = size(Ain);
    for icat=1:ncc
      dmm         = squeeze(Ain(icat,1:2:m4,1:2:n4));
      dmm(Ilo)    = 0;
      dmm(Ilnd)   = 0;
      A(icat,:,:) = dmm;
    end
  end

  nc_varput(fout,nname,A);

end

fprintf('%s is created\n',fout);
  




