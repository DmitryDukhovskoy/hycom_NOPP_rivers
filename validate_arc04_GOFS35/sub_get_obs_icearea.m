% Ice area for tot arctic and arctic ocean from bootstrap NSIDC NOAA
% see read_seaice_NSIDCdata.m
function [AA,TM] = sub_get_obs_icearea(d1,d2,fld) 

if isempty(fld)
	fld = 'TotArc';
end

pthiconc   = '/nexsan/people/ddmitry/Net_ocean/SeaIce/';
fmat = sprintf('%sBootstrap_Arctic_IceArea.mat',pthiconc);
load(fmat);
YR = ICE.Year;
MO = ICE.Month;
TotArc = ICE.TotArctArea_km2;
ArcOc  = ICE.ArctOcnArea_km2;

dv1 = datevec(d1);
yr1 = dv1(1);
mo1 = dv1(2);
dv2 = datevec(d2);
yr2 = dv2(1);
mo2 = dv2(2);

i1 = find(YR==yr1 & MO==mo1);
i2 = find(YR==yr2 & MO==mo2);
if isempty(i1) | isempty(i2),
	fprintf('sub_get_obs_icearea ERR: Could not locate dates\n');
	error('Stop')
end

if strncmp(fld,'TotArc',6)
	AA = TotArc(i1:i2);
else
	AA = ArcOc(i1:i2);
end

AA=AA(:);

icc=0;
for ii=i1:i2
	icc=icc+1;
	yr = YR(ii);
	mo = MO(ii);
	TM(icc,1) = datenum(yr,mo,1);
	TM(icc,2) = datenum(yr,mo,30);
end

return	

