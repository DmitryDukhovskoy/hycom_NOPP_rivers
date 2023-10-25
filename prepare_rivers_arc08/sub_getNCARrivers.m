function RVR = sub_getNCARrivers(RVR,Year,Mo);
% NCAR/UCAR rivers, monthly fields
% 1992 - 2014, then use 2014
pthncar = '/Net/ocean/ddmitry/UCAR_RIVERS/';
fnm = sprintf('%scoastal-stns-Vol-monthly.updated-Aug2014.nc',...
	      pthncar);
Time = nc_varget(fnm,'time');
YRS  = floor(Time./100);
MM   = Time - YRS*100;
DD   = MM*0+1;
TM   = datenum(YRS,MM,DD);

dnmb = datenum(Year,Mo,1);
ix = find(TM==dnmb);

% Find Arctic rivers, if needed
%if isempty(RVR) | ~isfield(RVR,'River_index')
ocn = nc_varget(fnm,'ocn_name');
[a1,a2] = size(ocn);
rivname = nc_varget(fnm,'riv_name');
cc = 0;
rindx = [];

for jk=1:a1
  ss = ocn(jk,1:3);
  if strncmp(ss,'ARC',3)
    cc = cc+1;
    rindx(cc,1)=jk;
  end
end

%end

RVR.Nrivers  = length(rindx);
RVR.Riv_name = rivname(rindx,1:15);

it1 = find(TM == datenum(1993,1,1));
it2 = find(TM == datenum(2014,12,1));
A = Flow(it1:it2,rindx);

Flow = nc_varget(fnm,'FLOW');
dmm = Flow(ix,rindx);  

return