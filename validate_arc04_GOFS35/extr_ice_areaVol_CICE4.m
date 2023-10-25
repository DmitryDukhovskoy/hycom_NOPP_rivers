% Extract time series of sea ice area & volume
% from specified experiments
% Old simulation with HYCOM-CICE4 GOFS3.1
% Used for tracer studies

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

pthout   = '/Net/kronos/ddmitry/hycom/ARCc0.04/datamat/cice_mnth/';
pthmat   = pthout;
pthtopo  = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthtopo4 = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';


f_add = 0;  % what experiment to add to existing saved, = 0 - all experiments

expt = 112;
expt1 = 'BL99Tfrz';


ftopo = sprintf('%sdepth_ARCc0.08_11.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
Lmsk = HH*0;
Lmsk(HH<0)=1;

[DX,DY]=sub_dx_dy(LON,LAT);
Acell08=DX.*DY; % Grid cell area, m2


% Output:
fmatout = sprintf('%s008hycom_cice4_112_seaice_VolArea.mat',pthmat);

if f_add >0; % add to existing saved experiments
  fprintf('Loading %s\n',fmatout);
  load(fmatout);
end

dd = 5; 
YR = 2017;
d1 = datenum(YR,4,1);
d2 = datenum(YR,12,31);

pthbin = sprintf('/nexsan/hycom/ARCc0.08_%3.3i/data/%4.4i_cice/',expt,YR);
if expt==110
  pthbin = sprintf('/nexsan/archive/ARCc0.08_110/data/%4.4i_cice/',YR);
end

Acell = Acell08;

nme = expt1;
ICE.Name = nme;

icc = 0;
for dnmb = d1:dd:d2
	dv   = datevec(dnmb);
	mo   = dv(2);
	mday = dv(3);

	fprintf('Reading %s %s\n',nme,datestr(dnmb));
	flnm = sprintf('%s%3.3i_cice_inst.%i-%2.2i-%2.2i-00000.nc',pthbin,expt,YR,mo,mday);

	if ~exist(flnm,'file');
		fprintf('Not found %s\n',flnm);
		continue;
	end

	fprintf('Reading %s \n',flnm);
	Hi = squeeze(nc_varget(flnm,'hi'));
	Ci = squeeze(nc_varget(flnm,'aice'));

	cice0=0.001;  % for ice extent - should be 0.15, for ice area - all ice >0% is considered
%    Ic = find(Ci>=0.15); % do not need this for the area, only for extent!
	Ic = find(Ci>=cice0);
	cci = Ci;
	cci(Ci<cice0)=0;
	Aice = nansum(Acell(Ic).*cci(Ic)*1e-6); % km2
	Vice = nansum(Acell(Ic).*cci(Ic).*Hi(Ic)*1e-9);  % km3

	fprintf('  A=%4.1fe6 km2, Vol=%6.4d km3\n',Aice*1e-6',Vice);
	icc = icc+1;
	ICE.TM(icc)       = dnmb;
	ICE.Area_km2(icc) = Aice;
	ICE.Vol_km3(icc)  = Vice;

%keyboard

end

fprintf('Saving %s\n',fmatout);
save(fmatout,'ICE');


