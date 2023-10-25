function [TM, Rgr] = sub_read_Greenland_v3;
% Total greenland runoff for specified year, month
% and cumulative FWF anomaly up to this date
%
% Read Greenland runoff - Bamber data - updated version of 2018
% with 2011-2016 
% calculate surplus flux
% wrt pre-1990 mean
%dv0 = datevec(dnmb);
%dnmb = datenum(dv0(1),dv0(2),1);
fprintf('Reading Greenland runoff data ...\n');

PTH.river   = '/nexsan/people/ddmitry/Net_ocean/arctic_AOregimes/data/GreenlandRunoffv3/';
% Read Bamber's data: fluxes = km3/mo
% Total FW flux = D+Rg+Rt (no CAA, Svlabard here);
% Shown in Bamber's Fig.3, 2018
fnm = sprintf('%sFWF17.v3.nc',PTH.river);
tmm = double(nc_varget(fnm,'TIME'));
Xgr = nc_varget(fnm,'lon');
Ygr = nc_varget(fnm,'lat');
Rt  = double(nc_varget(fnm,'runoff_tundra')); % tundra runoff, km3/mo
Rg  = double(nc_varget(fnm,'runoff_ice')); % GrIS runoff - meltwater
Ds  = nc_varget(fnm,'solid_ice'); % solid ice
LGr = nc_varget(fnm,'LSMGr'); % Greenland mask
TM  = datenum(1958,1,1)+tmm; 
DV  = datevec(TM);
nrc = length(DV);


Rgr = [];
for it=1:nrc
  if mod(it,20)==0,
    fprintf('  %6.2f%% ...\n',it/nrc*100);
  end
  
  dmm  = squeeze(Ds(it,:,:));
  Di   = abs(dmm).*LGr; % Solid disch, km3/mo
  dmm  = squeeze(Rt(it,:,:));
  Rti  = abs(dmm).*LGr; % Gr tundra
  dmm  = squeeze(Rg(it,:,:));
  Rgi  = abs(dmm).*LGr; % Greenland meltwater

  FWF  = Di+Rti+Rgi; % total FWF =  disch+tundra+GrMelt
  %    Rgr=GR(k).runoff;
%  Rgr(it,1)  = nansum(nansum(FWF))*1e9/(3600*24*mday(k));  % km3/mo -> m3/s  
  Rgr(it,1)  = nansum(nansum(FWF));  % km3/mo  
end


return
