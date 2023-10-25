% Correct concentration of Bering tracers
% extracted in extr_trc_day.m
%
% Saved monthly mean tracer concentrations
% integrated over some depth:
% 0 - 50 m - mixed layer
% 150 - 50 m
%
% river passive tracers 
% distributed along the Greenland Coast
%
% Adjust Bering Strait tracer concenctration:
% Distribute FW flux (0.08-0.1 Sv) over all points where Tracer=1
% 0.08 Sv = 2522.8 km3/yr
% 0.1 Sv = 3200 km3/yr
% Model long-term mean is 2778.9 km3/yr (anls_fluxes/calc_FWFlux_Bering.m)
%
%
% 

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers
startup;

close all
clear

s_mat = 1;

cBr = sub_BeringTrConc(0);% coefficient for Bering Str., tracer, kg/m3

rg = 9806; 
%cBr = 14; % coefficient for Bering Str. tracer - not right
%cBr0 = 0.1; % old cft =1 if not changed or whater it was assigned 
            % in extr_tracer_* 
%crct = cBr/cBr0; % correction
nTr = 5; 

%yr0 = 2005; 

lr1 = [0,-50];
lr2 = [-50,-150];

yr1=1993;
yr2=2015;

fprintf('Years to extract: %i-%i\n',yr1,yr2);

% Experiments:
% 102 - test simulation, HYCOM GLBb0.08 nest: 1 file kept constant (1993,1,1)

regn = 'ARCc0.08';
expt = 110;  
%pthfig  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/fig_trac/',expt);
pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);


ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

hmsk=HH;
hmsk(HH<0)=nan;

xlim1 = 0;
xlim2 = nn;
ylim1 = 0;
ylim2 = mm;

% Read fields:
ip1=1;
mold = 0;
dday = 7; 
for iyr = yr1:yr2
  for imo = 1:12
    fprintf('Changing tracer %i, %i/%2.2i\n',nTr,iyr,imo);
    
    fmat = sprintf('%strcr_dpthav_%i%2.2i.mat',pthmat,iyr,imo);
    fprintf('Loading %s\n',fmat);
    load(fmat);
    
    for ilv=1:2
      Tr = squeeze(TRCR(ilv).TR(nTr,:,:));
      tmx = max(max(Tr)); % max conc
      crct = cBr/tmx; % make max conc = cBr
%      keyboard
      Tr = Tr*crct;
      TRCR(ilv).TR(nTr,:,:)=Tr;
    end
    
   if s_mat>0
     fprintf('Rewriting %s\n',fmat);
     save(fmat,'TRCR','-v7.3');
   end;
    
  end % month  
end  % year


%keyboard



