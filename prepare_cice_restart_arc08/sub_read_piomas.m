function AP = sub_read_piomas(fld,yr,imo);
% Read PIOMAS sea ice fields
% For daily data imo = month*100 + day
iday=0;
if imo<13
  imo=imo;
else
  dmm=imo;
  imo=floor(dmm/100);
  iday=dmm-imo*100;
end

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
switch(fld);
 case('ithkn');
  fthk=sprintf('%sheff.H%i',PTH.data,yr);
  fid=fopen(fthk,'r','ieee-le');
% Read 1st month:
  frewind(fid);
  for im=1:imo
    aa=fread(fid,nxy,'float32');
    ithk=reshape(aa,[nxl,nyl])';
  end;
  AP.Fld   = fld;
  AP.File  = fthk;
  AP.Field = ithk;
  AP.LON   = LONs;
  AP.LAT   = LATs;

 case('hiday'); % daily thickness
  fthk=sprintf('%shiday.H%i',PTH.data,yr);
  fid=fopen(fthk,'r','ieee-le');
% Read 1st month:
  frewind(fid);
  for im=1:imo
    for idd=1:iday
      aa=fread(fid,nxy,'float32');
      ithk=reshape(aa,[nxl,nyl])';
    end
  end;
  AP.Fld   = fld;
  AP.File  = fthk;
  AP.Field = ithk;
  AP.LON   = LONs;
  AP.LAT   = LATs;



 case('itemp');
  fthk=sprintf('%stice0.H%i',PTH.data,yr);
  fid=fopen(fthk,'r','ieee-le');
% Read 1st month:
  frewind(fid);
  for im=1:imo
    aa=fread(fid,nxy,'float32');
    ithk=reshape(aa,[nxl,nyl])';
  end;
  AP.Fld   = fld;
  AP.File  = fthk;
  AP.Field = ithk;
  AP.LON   = LONs;
  AP.LAT   = LATs;

 case('snow');
  fthk=sprintf('%ssnow.H%i',PTH.data,yr);
  fid=fopen(fthk,'r','ieee-le');
% Read 1st month:
  frewind(fid);
  for im=1:imo
    aa=fread(fid,nxy,'float32');
    ithk=reshape(aa,[nxl,nyl])';
  end;
  AP.Fld   = fld;
  AP.File  = fthk;
  AP.Field = ithk;
  AP.LON   = LONs;
  AP.LAT   = LATs;
 
 case('area');
% Read ice conc
  fthk=sprintf('%sarea.H%i',PTH.data,yr);
  fid=fopen(fthk,'r','ieee-le');
% Read 1st month:
  frewind(fid);
  for im=1:imo
    aa=fread(fid,nxy,'float32');
    ithk=reshape(aa,[nxl,nyl])';
  end;
  AP.Fld   = fld;
  AP.File  = fthk;
  AP.Field = ithk;
  AP.LON   = LONs;
  AP.LAT   = LATs;
 
 case('aiday');  % Daily ice area
% Read ice conc
  fthk=sprintf('%saiday.H%i',PTH.data,yr);
  fid=fopen(fthk,'r','ieee-le');
% Read 1st month:
  frewind(fid);
  for im=1:imo
    for idd=1:iday
      aa=fread(fid,nxy,'float32');
      ithk=reshape(aa,[nxl,nyl])';
    end
  end;
  AP.Fld   = fld;
  AP.File  = fthk;
  AP.Field = ithk;
  AP.LON   = LONs;
  AP.LAT   = LATs;
 
  
end

return

