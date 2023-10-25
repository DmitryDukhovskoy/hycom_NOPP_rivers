% Plot Tracer fraction in the BG from 
% monthly-mean, depth-integrated 
% tracer concentrations 
% The Data are extracted in extr_trcr_mnth.m
% Option: plot anomalies relative to reference Year/month
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers
startup;

close all
clear

crctB = 0.296; % temporarily, correct too high estimate of Bering tracer
s_fig = 0;
s_mat = 2; % =0 - do not save mat file
           % = 1 donwl. tracer conc 0-50, 150-50 levels, extract for region, save
	   % = 2 skip extraction, only plot output
%dref = datenum(2006,1,15); % reference year/month, = 0 - no reference 
%ilev = 1; %=1 - mixed layer, 0-50; =2 - 50-150m
yr1 = 1993;
yr2 = 2016;

regn = 'ARCc0.08';
expt = 110;  
pthfig  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/fig_BGbudget/',expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
fmout   = sprintf('%sfraction_trcr_BG.mat',pthmat);

ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

% Grid cell spacing
[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2


xlim1 = 20;
xlim2 = nn;
ylim1 = 5;
ylim2 = 2000;

hmsk=HH;
hmsk(HH<0)=nan;

% Define BG:
iBG = [483, 1532;...
       684, 1345; ...
       757, 1423;...
       601, 1637;...
       483, 1532];


[XX,YY] = meshgrid((1:nn),(1:mm));
INP = inpolygon(XX,YY,iBG(:,1),iBG(:,2));
IN = find(INP==1);

chck=0;
if chck>0
  figure(10); clf;
  contour(HH,[0 0],'k');
  hold on;
  contour(HH,[-5000:1000:-10],'b');
  plot(iBG(:,1),iBG(:,2),'r.-');
  contour(LON,[-153 -153],'g');
  contour(LON,[-130 -130],'g');
  contour(LAT,[73 73],'g');
  contour(LAT,[83 83],'g');
  [JJ,II]=ind2sub(size(HH),IN);
  plot(II,JJ,'y.');
  axis('equal');
end

nint = 360;
c1 = -1;
c2 = 1;
CMP = create_colormap2_3(nint,c1,c2);
cmp = CMP.colormap;
cmp = smooth_colormap(cmp,18);
cmp(1,:) = [1 1 1];
cnt = CMP.intervals;

if s_mat<2
  cc=0;
  clear sumT
%  yr1 = 1993;
%  yr2 = 2016;
  for iyr=yr1:yr2
    for imo=1:12
      fmat = sprintf('%strcr_dpthav_%i%2.2i.mat',pthmat,iyr,imo);
      fprintf('Loading %s\n',fmat);
      if exist(fmat,'file')
	load(fmat);
      else
	fprintf(' =========  MISSING %s\n',fmat);
	cc=cc+1;
	for nTr=1:5
	  for ilv=1:2
	    sumT(cc,ilv,nTr) = nan;
	  end
	end
	continue 

      end



      cc=cc+1;
      for nTr = 1:5
	rr = [];
	for ilv=1:2
	  fprintf('Reading: %i/%2.2i, Tracer %i, Lev %i\n',...
		  iyr,imo,nTr,ilv);
	  Tr = squeeze(TRCR(ilv).TR(nTr,:,:));
	  Tr(Tr<=0)=nan;
  % Integrate Tracer over 1 region of interest per layers
  % To get average Tracer concentration within the region
	  dz = abs(TRCR(ilv).depth_av(2)-TRCR(ilv).depth_av(1));
          volR = sum(Acell(IN)*dz); % total vol, region
	  mnt = nansum(Tr(IN).*Acell(IN)*dz)/volR; % avrg Tr conc in region
	  mnTR(cc,ilv,nTr) = mnt; % avrg Tr conc in region
	  
  % Integrate Tracer over the whole domain and all depths - total mass
  % rr = mass tracer in the domain / total mass of the tracer
	  intgrTr_dom = nansum(nansum(Tr.*Acell*dz));
	  intgrTr_rgn = nansum(nansum(Tr(IN).*Acell(IN)*dz));
	  rr=intgrTr_rgn/intgrTr_dom; 
	  sumTdom(cc,ilv,nTr) = rr;
 %if iyr>1995; keyboard; end
	end

      end

    end % month
  end   % year
  if s_mat==1
    fprintf('Saving %s\n',fmout);
    save(fmout,'mnTR','sumTdom');
  end
else
  fprintf('Loading %s\n',fmout);
  load(fmout);
  
  mnTR(:,:,5)=mnTR(:,:,5)*crctB;
  sumTdom(:,:,5)=sumTdom(:,:,5)*crctB;
  
end
  


%keyboard

% sum over all tracers:
clear atot
for ilv=1:2
  dmm = squeeze(mnTR(:,ilv,:));
  smm = sum(dmm,2);
  atot(:,ilv) = smm;
end

% Ratio: Tracer mass/all tracers inside the region
clear R
for ilv=1:2
  dmm = squeeze(mnTR(:,ilv,:));
  for nTr=1:5
    amm = dmm(:,nTr);
    rr = amm./atot(:,ilv);
    R(:,ilv,nTr) = rr;
  end
end

nrc = size(R,1);
yrs=[0:nrc-1]/12+yr1;
%
% ----------------------
%

% ==========================================================
%   Plot fraction of tracers in the BG wrt to the
% total tracer mass in the domain
% by depth intervals (saved in the extracted time ser.)
% ==========================================================
figure(1); clf;
ilv=1;
dmm = squeeze(sumTdom(:,ilv,2));
dmm2 = squeeze(sumTdom(:,ilv,5));

% Mackenzie River
axes('position',[0.08 0.6 0.85 0.35]);
plot(yrs,dmm);
hold
plot(yrs,dmm2,'r');
yl1=0;
yl2=1.15*max(max(dmm));
set(gca,'tickdir','out',...
	'ylim',[yl1 yl2],...
	'xlim',[2003 max(yrs)],...
	'xgrid','on',...
	'xtick',[yr1:yr2],...
	'xminortick','on',...
	'yminortick','on',...
	'ygrid','on');
title('Mackenzie (b) & Bering W, Fraction in the BG, HYCOM, Mixed L.');

% Eurasian Rivers
axes('position',[0.08 0.1 0.85 0.35]);
ilv=1;
it=2;
dmm = squeeze(sumTdom(:,ilv,3));
dmm2 = squeeze(sumTdom(:,ilv,4));
plot(yrs,dmm); % E. Eurasian
hold
plot(yrs,dmm2,'r'); % W. Eurasian
yl1=0;
yl2=1.15*max(max(dmm));
set(gca,'tickdir','out',...
	'ylim',[yl1 yl2],...
	'xlim',[2003 max(yrs)],...
	'xgrid','on',...
	'xtick',[yr1:yr2],...
	'xminortick','on',...
	'yminortick','on',...
	'ygrid','on');
title('E.(b) and W.EurasianR., Fraction in the BG, HYCOM, Mixed L.');

btxt='fraction_trcr_BG.m';
bottom_text(btxt,'pwd',1);



% Fraction of Tracer: - Mean concentration in the BG 
figure(2); clf;

% Plot individual tracers:
axes('position',[0.08 0.8 0.85 0.15]);
dmm = squeeze(mnTR(:,ilv,2)); % Mackenzie
%dmm2 = squeeze(sumT(:,ilv,5)); % Bering Str. Water
%dmm = squeeze(R(:,ilv,2))*100; % Mackenzie
%dmm2 = squeeze(R(:,ilv,5))*100; % Bering Str. Water
plot(yrs,dmm,'linewidth',2);
%hold 
%plot(yrs,dmm2,'r');

yl1=0.95*(min(dmm));
yl2=1.1*max(dmm);
set(gca,'tickdir','out',...
	'ylim',[yl1 yl2],...
	'xlim',[2003 max(yrs)],...
	'xtick',[yr1:yr2],...
	'xgrid','on',...
	'ygrid','on',...
	'Fontsize',14,...
	'xminortick','on',...
	'yminortick','on');
title('Mackenzie, kg/m3');

% 
axes('position',[0.08 0.55 0.85 0.15]);
dmm = squeeze(mnTR(:,ilv,3)); % E. Eurasian Rivers
plot(yrs,dmm,'linewidth',2);
yl1=0.95*(min(dmm));
yl2=1.1*max(dmm);
set(gca,'tickdir','out',...
	'ylim',[yl1 yl2],...
	'xlim',[2003 max(yrs)],...
	'xtick',[yr1:yr2],...
	'xgrid','on',...
	'ygrid','on',...
	'Fontsize',14,...
	'xminortick','on',...
	'yminortick','on');
title('E.Eurasian Rivers, kg/m3');

axes('position',[0.08 0.33 0.85 0.15]);
dmm = squeeze(mnTR(:,ilv,4)); % W. Eurasian Rivers
plot(yrs,dmm,'linewidth',2);
yl1=0.95*(min(dmm));
yl2=1.1*max(dmm);
set(gca,'tickdir','out',...
	'ylim',[yl1 yl2],...
	'xlim',[2003 max(yrs)],...
	'xtick',[yr1:yr2],...
	'ygrid','on',...
	'xgrid','on',...
	'xminortick','on',...
	'yminortick','on');
title('W.Eurasian Rivers,kg/m3');

axes('position',[0.08 0.1 0.85 0.15]);
dmm = squeeze(mnTR(:,ilv,5)); % Pacific
plot(yrs,dmm,'linewidth',2);
yl1=0.95*(min(dmm));
yl2=1.1*max(dmm);
set(gca,'tickdir','out',...
	'ylim',[yl1 yl2],...
	'xlim',[2003 max(yrs)],...
	'xtick',[yr1:yr2],...
	'xgrid','on',...
	'ygrid','on',...
	'Fontsize',14,...
	'xminortick','on',...
	'yminortick','on');
title('Pacific Water, kg/m3');



btxt='fraction_trcr_BG.m';
bottom_text(btxt,'pwd',1);

if s_fig >0
  fgnm=sprintf('%sarc008_110_fraction_trcrBG',pthfig);
  fprintf('Saving %s\n',fgnm);
  print('-dpng','-r300',fgnm);
end

% Plot Mackenzie and Pacific Water, kg/m3
figure(4); clf
axes('position',[0.08 0.55 0.85 0.35]);
dmm = squeeze(mnTR(:,ilv,2)); % Mackenzie
dmm2 = squeeze(mnTR(:,ilv,5)); % Bering Str. Water
%dmm = squeeze(R(:,ilv,2))*100; % Mackenzie
%dmm2 = squeeze(R(:,ilv,5))*100; % Bering Str. Water
plot(yrs,dmm,'linewidth',2);
hold 
plot(yrs,dmm2,'r','linewidth',2);

yl1=0.95*(min([dmm;dmm2]));
yl2=1.1*max([dmm;dmm2]);
set(gca,'tickdir','out',...
	'ylim',[yl1 yl2],...
	'xlim',[2003 max(yrs)],...
	'xtick',[yr1:yr2],...
	'xgrid','on',...
	'ygrid','on',...
	'xminortick','on',...
	'yminortick','on',...
	'Fontsize',14);
title('Mackenzie (b) and Pacifc W.(r), kg/m3');

btxt='fraction_trcr_BG.m';
bottom_text(btxt,'pwd',1);

