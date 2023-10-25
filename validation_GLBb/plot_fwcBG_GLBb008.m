% Time series extracted from
% HYCOM GLBb0.08
% see: hycom_NOPP_rivers/validation_GLBb/fwcBG_GLBb008.m
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

s_fig  = 0;
fld    = 'salt';

fprintf('saving mat = %i\n',s_mat);

rg=9806;  % convert pressure to depth, m
Sref=34.8; % N.Atl. is too saline
TV = '07';

pthARC = '/Net/mars/ddmitry/hycom/GLBb2ARCc0.08/mean/';
pthmat = '/Net/mars/ddmitry/hycom/GLBb0.08/data_mat/';
pthfig = '/Net/mars/ddmitry/hycom/GLBb0.08/fig_fw/';
%ftopo = sprintf('%sGLBb_T07_subset_Natl.mat',pthmat);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
%pthglb = '/nexsan/GLBb0.08/GLBb0.08_191/topo/';

TM = [];
Fwc = [];
for yr = 2009:2017
  fmat = sprintf('%sFWC_BGvol_GLBb008_%i.mat',pthmat,yr);
  fprintf('Loading %s\n',fmat);
  load(fmat);
  tm = FWC.TM;
  tm = tm(:);
  TM = [TM;tm];
  fwc = FWC.Fwc_km3*1e-3; % thousand km3
  fwc = fwc(:);
  Fwc = [Fwc;fwc];
end

nrc = length(TM);
npyr= round(365/5); % rec/year
DV  = datevec(TM);
yr1 = DV(1,1);
yr2 = DV(end,1);
yrs = [0:nrc-1]/npyr+yr1;

cc=0;
for iyr=2009:0.25:2018
  cc=cc+1;
  if round(iyr)~=iyr
    xtcklb{cc}=' ';
  else
    xtcklb{cc}=sprintf('%i',round(iyr));
  end
end

% Delete suspecious jumps:
Fwc(27:45)=nan;
Fwc(367) = nan;
Fwc(441) = nan;
Fwc(515) = nan;

figure(1); clf;
axes('Position',[0.08 0.45 0.85 0.45]);
plot(yrs,Fwc);
set(gca,'tickdir','out',...
	'ylim',[15 21],...
	'xlim',[2009 2017.5],...
	'xtick',[2009:0.25:2018],...
	'xticklabel',xtcklb);
title('FWC, x1e3 km^3, Beauf.Gyre, GLBb0.08 HYCOM+NCODA GOFS3.1',...
      'Fontsize',11);

btxt = 'plot_fwcBG_GLBb0.08';
bottom_text(btxt,'pwd',1,'position',[0.02 0.3 0.7 0.1]);

if s_fig>0
  fgnm = sprintf('%sBG_FWC_GLBb008_%i-%i',pthfig,yr1,yr2);
  fprintf('Saving %s\n',fgnm);
  print('-dpng','-r250',fgnm);
end




