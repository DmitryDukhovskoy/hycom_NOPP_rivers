% Read PIOMAS sea ice fields
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close

PTH.topo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
PTH.data = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/PIOMAS_ice_data/';
PTH.rest = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/data_ice_restart/';

fgrds = sprintf('%sgrid.dat',PTH.data); % grid for scalar fields
fgrdv = sprintf('%sgrid.dat.pop',PTH.data); % grid for vect. fields

nxl=360;
nyl=120;
nxy=nxl*nyl;

% read lon/lat scalar fields
dmm = load(fgrds);
[a1,a2]=size(dmm);
nrw=nxy/a2;

LONs = dmm(1:nrw,:);
LONs = reshape(LONs',[nxy,1]);
LONs = reshape(LONs,[nxl,nyl])';

LATs = dmm(nrw+1:end,:);
LATs = reshape(LATs',[nxy,1]);
LATs = reshape(LATs,[nxl,nyl])';

% read lon/lat vector fields
dmm = load(fgrdv);
[a1,a2]=size(dmm);
nrw=nxy/a2;

LONv = dmm(1:nrw,:);
LONv = reshape(LONv',[nxy,1]);
LONv = reshape(LONv,[nxl,nyl])';

LATv = dmm(nrw+1:nrw+nrw,:);
LATv = reshape(LATv',[nxy,1]);
LATv = reshape(LATv,[nxl,nyl])';

% Read sea ice thickness
% only monthly mean
fthk=sprintf('%sheff.H1993',PTH.data);
fid=fopen(fthk,'r','ieee-le');
% Read 1st month:
frewind(fid);
for im=1:1
  aa=fread(fid,nxy,'float32');
  ithk=reshape(aa,[nxl,nyl])';
end;



