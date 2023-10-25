% Time series of tracer mass by layers
% The Tracer Data are extracted in extr_MassTrcr_mnth.m
% mat files: MassTr01_lrs_*.mat for Greenland tracer

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers
startup;

close all
clear

nTr = 1; 
yr1 = 1993;
yr2 = 2016;
nbx = 5; % plot boxes =1,..., nbx

s_fig  = 0;

% Specify levels to plot
ilv = 1; % 0-50m
%ilv = 2; % 50-150m
%ilv = 3; % 150-300m
%ilv = 4; % 300-500 m
%ilv = 5; % whole depth 

% Vertical layers
LRS = load('LRS.dat');
nlrs= length(LRS); 


zz1 = LRS(ilv,1);
zz2 = LRS(ilv,2);

dz=abs(zz2-zz1);

fprintf('Time Ser. Tracer Mass budget, ilv=%i, %i - %i, nTr=%i\n',ilv,zz1,zz2,nTr);


regn = 'ARCc0.08';
expt = 110;  
pthfig  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/fig_NAtl/',expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
pthsav  = sprintf('/Net/tholia/ddmitry/hycom/ARCc0.08/data_mat%3.3i/',expt);
%fmout   = sprintf('%sGrTrFrct_GrWVol_Regions_NAtl_lev%2.2i.mat',pthmat,ilv);
btx     = 'get_tser_trmass.m';

LRS = load('LRS.dat');
nlrs= length(LRS); 

ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
IN  = find(HH<0);

% Grid cell spacing
[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2


cc=0;
TRM = [];
%  clear sumT
for iyr=yr1:yr2
  for imo=1:12
%pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
% fmat = sprintf('%sMassTr%2.2i_lrs_%i%2.2i.mat',pthmat,nTr,iyr,imo);
    dnmb = datenum(iyr,imo,15);
    dv0  = datevec(dnmb);
    iyr  = dv0(1);
    iday = dnmb-datenum(iyr,1,1)+1;
    imo  = dv0(2);

    fmat = sprintf('%sMassTr%2.2i_lrs_%i%2.2i.mat',pthmat,nTr,iyr,imo);
    fprintf('Loading %s\n',fmat);
    if exist(fmat,'file')
      load(fmat);
    else
      fprintf(' =========  MISSING %s\n',fmat);
      return
    end
    fprintf('Reading: %i/%2.2i, Tracer %i, Lev %i\n',...
	    iyr,imo,nTr,ilv);
    cc=cc+1;
    for ilv=1:5
      Tr = squeeze(TRCR(ilv).MassTr_kg);
      Tr(Tr<1e-23)=nan;

      TRM(ilv,cc)=nansum(nansum(Tr));
    end
    
  end
end


YRS = [1993:1/12:yr2+0.99];
figure(1); clf;
hold
for ilv=5:5
  trm = TRM(ilv,:);
  plot(YRS,trm);
end

set(gca,'xlim',[yr1 yr2+1],...
	'xgrid','on',...
	'ygrid','on');

% Compare with cumulative Gr FWF anomaly
frv = sprintf('%sGreenland_cumFWFlux.mat',pthmat);
fprintf('Loading %s\n',frv);
load(frv);  % cFWF, km3



% Whole depth
trm = TRM(5,:);
dmm = reshape(trm,[12,24]);
dmm = mean(dmm);
RT = dmm/max(dmm)*max(cFWF);

nyr = 24;
YY = [1993:2016];

% Estimated FW flux in the experiment
% when a constant GrFWF anomaly is 
% imposed: 211.7 km3/yr
cFa = [1:24]*208.6;

figure(2); clf;
axes('Position',[0.08 0.45 0.85 0.45]);
hold
hb = bar(Ygr,cFWF);
set(hb,'Facecolor',[0.8 0.8 0.8]);
plot(YY,RT,'k','Linewidth',2);
plot(YY,cFa,'k--','Linewidth',2);


title('ARCc0.08-110, Vol (km3) of GrFW anom from total TrMass and GrRunoff');


set(gca,'tickdir','out',...
        'xlim',[yr1-0.2 yr2+0.5],...
	'ylim',[-200 5300],...
	'ytick',[-500:500:6000],...
	'xtick',[1993:2017],...
	'xgrid','on',...
	'ygrid','on');

bottom_text(btx,'pwd',1,'Position',[0.02 0.3 0.4 0.1]);

