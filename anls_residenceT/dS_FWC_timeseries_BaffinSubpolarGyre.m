% Time series of FWC and dltS 
% in Subpolar Gyre region
%
% Uses corrected code extracting tracer
% that integrates tracer over the specified layers
% see: extr_MassTrcr_month
% Estimate S change in the upper layer due to surplus FW flux
% from Greenland 
% Note estimate of Greenland runoff fraction
% for specified boxes may use wrong
% code - tracer concentration averaged over layers
% should be integrated giving a total mass within
% a layer

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers
startup;

close all
clear

s_fig = 0;
f_extr=0;

nTr = 1; 
IntSrf = 0; % = 1 - integrates over all layers from 0 - ilv, for ilv=1 is the same
%ilv = 1; % 0-50m
%ilv = 2; % 50-150m
%ilv = 3; % 150-300m
%ilv = 4; % 300-500 m
ilv = 5; % whole depth <- do not use this

% Vertical layers
LRS = load('LRS.dat');
nlrs= length(LRS); 


zz1 = LRS(ilv,1);
zz2 = LRS(ilv,2);

if IntSrf==1; zz1=0; end;

dz=abs(zz2-zz1);

fprintf('FW Volume Subpolar Gyre, ilv=%i, %i - %i, nTr=%i, 1993-2016\n',ilv,zz1,zz2,nTr);

regn = 'ARCc0.08';
expt = 110;  
pthfig  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/fig_NAtl/',expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
pthsav  = sprintf('/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_mat%3.3i/',expt);
%fmout   = sprintf('%sfraction_trcr_NAtl.mat',pthmat);

ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

% Grid cell spacing
[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2
[XX,YY] = meshgrid((1:nn),(1:mm));

hmsk=HH;
hmsk(HH<0)=nan;

%
hZ = abs(LRS(ilv,2)-LRS(ilv,1));

if IntSrf==1
  hZ = dz;
end


%dSregions_map_v2


% Subpolar Gyre:
% Labr - Irm - eastern & central N. Atlantic
% Use Davis Str and Denmark Str - same locations
% as in the POP gates to calculate
% fluxes 
  IGR = [  430  729
	 427  957
	 459  1009
	 556  1066
	 739  1092
	 748  908
         810  621
         855  530
         930         501
        1031         445
        1099         410
        1192         403
        1228         321
        1218         185
        1080         145
         444         145
         416         246
         367         492
         374         591
         379         729];

% case('spg2')
%
%end



INP = inpolygon(XX,YY,IGR(:,1),IGR(:,2));
INGr= find(INP==1 & HH<0);
II = find(HH<0);



frv = sprintf('%sGreenland_cumFWFlux.mat',pthmat);

f_griv = 0; %=1 - rederive cumulative Greenland FWFlux from Bamber
if f_griv==1
  [TMg,Fgr]=sub_read_Greenland_v3; % Greenland total FWF, km3/mo
  % Calculate annual anomalies
  % and cumulative up to date
  DVg = datevec(TMg);
  ngr = length(TMg);
  nyr =ngr/12;
  dmm = reshape(Fgr,[12,nyr]);

  Fyr = sum(dmm); % km3/yr
  Ygr = [DVg(1,1):DVg(end,1)];
  ii=find(Ygr==1990);
  Fmn = mean(Fyr(1:ii));
  ism = find(Ygr==dv0(1,1));
  cFWF = cumsum(Fyr-Fmn);
  fwf0 = cFWF(ism); % km3
  save(frv,'cFWF','Ygr');
else
  fprintf('f_griv %i, Loading %s\n',f_griv,frv);
  load(frv);
%  ii=find(Ygr==1990);
end  


% Tracer fraction in grid cells
fextr=sprintf('%sFWC_BaffSbpGyre_TimeSer.mat',pthmat);
if f_extr==1
  nTr = 1;
  iyy=0;
  for iyr=1993:2016
    iyy=iyy+1;
    for im=1:12
      dnmb = datenum(iyr,im,15);
      dv0  = datevec(dnmb);
      yr   = dv0(1);
      iday = dnmb-datenum(yr,1,1)+1;

      ism = find(Ygr==dv0(1,1));
      fwf0 = cFWF(ism); % km3

      rr = sub_fraction_tracerMass(expt,regn,Acell,HH,dnmb,nTr,ilv,hZ,IntSrf);

  % Estimate volume of Greenland surplus FW in grid cell
  % = ratio of FW tracer*total FWF (integrated over time 1990-date requested)
      Vfw = fwf0*rr*1e9; % m3 Gr FWF anomaly integrated in time
  % Integrate vol of FW in the region
      vfw=nansum(Vfw(INGr))*1e-9; % m3->km3
      VFW(im,iyy)=vfw;

    end
  end
  fprintf('Saving %s\',fextr);
  save(fextr,'VFW');
else
  fprintf('Loading %s\n',fextr);
  load(fextr);
end


[m,n]=size(VFW);
VFWy=mean(VFW);
VF=reshape(VFW,[m*n,1]);

YR=[1993:2016];
figure(1); clf;
axes('Position',[0.1 0.5 0.85 0.42]);
plot(YR,VFWy,'Linewidth',2);
stl=sprintf('Tracer-based Vol GrFW anomaly in Baffin+SbpGyre, max=%6.1f km^3',...
      max(VFWy));
title(stl);

set(gca,'tickdir','out',...
	'xlim',[1993 2016],...
	'ylim',[-100 2500],...
	'xgrid','on',...
	'ygrid','on',...
	'Fontsize',14);
btx='dS_FWC_timeseries_BaffinSubpolarGyre.m';
bottom_text(btx,'pwd',1,'Position',[0.08 0.4 0.4 0.05]);

axes('Position',[0.6 0.05 0.3 0.3]);
hold on
plot(XX(INGr),YY(INGr),'g.');
contour(HH,[0 0],'k');
axis('equal');
set(gca,'xlim',[300 1300],...
	'ylim',[50 1300],...
	'xtick',[],...
	'ytick',[]);

    
