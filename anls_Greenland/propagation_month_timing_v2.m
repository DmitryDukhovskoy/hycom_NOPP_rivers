% Estimate timing/propagation rate
% of the Greenland FW 
% using tracer
% Here: 
% Time of reaching 80 or 90% of max ratio = m(i,j)/Mtr overall
%
% use extracted time series, monthly  mean fields
% within the vertical layers
%
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

rmax = 0.8; % time when tr>rmax
s_fig  = 0;
f_maxc = 1; % derive max conc in grid cells
f_maxt = 1; % derive time of max conc in grid cells
%dref = datenum(2006,1,15); % reference year/month, = 0 - no reference 
%ilev = 1; %=1 - mixed layer, 0-50; =2 - 50-150m

% Specify levels:
%ilv = 1; % 0-50m
%ilv = 2; % 50-150m
%ilv = 3; % 150-300m
ilv = 4; % 300-500 m
%ilv = 5; % whole depth <- do not use this

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
fmout   = sprintf('%sGrTrFrct_GrWVol_Regions_NAtl_lev%2.2i.mat',pthmat,ilv);
btx     = 'propagation_month_timing_v2.m';

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


flcmx=sprintf('%smax_trmass_v2_mo_lev%2.2i.mat',pthmat,ilv);
if f_maxc==1
  cc=0;
  CMX = HH*0;
%  clear sumT
  for iyr=yr1:yr2
    cc=0;
    dmm=HH*0;
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
      Tr = squeeze(TRCR(ilv).MassTr_kg);
      MTr = nansum(nansum(squeeze(TRCR(5).MassTr_kg)));
      rr = Tr./MTr*1e12;
%      dmm= dmm+Tr;
%      Tr = dmm/cc;
      CMX = max(CMX,rr);
    end
  end
  
  CMX(HH>=0)=nan;
  CMX(CMX<1e-10)=nan;
  
%  flcmx=sprintf('%smax_trmass_lev%2.2i.mat',pthmat,ilv);
  fprintf('Saving %s\n',flcmx);
  save(flcmx,'CMX');
end

fprintf('Loading %s\n',flcmx);
load(flcmx);


% Get time 
TMX = HH*0;
tm0 = datenum(yr1,1,15);
fltmx=sprintf('%sTmax_trmass_v2_mo_lev%2.2i_r%3.3i.mat',pthmat,ilv,rmax*100);
if f_maxt==1
  for iyr=yr1:yr2
    cc=0;
    dmm=HH*0;
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
      Tr = squeeze(TRCR(ilv).MassTr_kg);
      MTr = nansum(nansum(squeeze(TRCR(5).MassTr_kg)));
      rr = Tr./MTr*1e12;
%      cc=cc+1;
%      dmm = dmm+Tr;
%      Tr = dmm/cc;
    
      II = find(rr>=rmax*CMX & TMX==0);
      if ~isempty(II),
        TMX(II)=dnmb-tm0;
      end
      fprintf('Finding time %i/%2.2i, # points found: %i\n',...
	      iyr,imo,length(II));
    end
    
  end
  TMX(HH>=0)=nan;
  TMX(TMX<=1)=nan;
  fprintf('Saving %s\n',fltmx);
  save(fltmx,'TMX');
end

fmat = sprintf('%sMassTr%2.2i_lrs_%i%2.2i.mat',pthmat,nTr,2016,7);
fprintf('Loading %s\n',fmat);
load(fmat);
fprintf('Reading: %i/%2.2i, Tracer %i, Lev %i\n',...
	2016,12,nTr,ilv);
Tr = squeeze(TRCR(ilv).MassTr_kg);
mtr=nanmean(nanmean(Tr));
I=find(Tr<0.1*mtr);


fprintf('Loading %s\n',fltmx);
load(fltmx);
TMX(I)=nan;
TMX=TMX./365.25;

% Spatial filtering:
sgmx = 2;
sgmy = sgmx;
npnts = 4*sgmx; % # of pnts in 1 dir
ilm1 = 10;
ilm2 = 1550;
jlm1 = 1;
jlm2 = 1600;
TMXf = sub_gauss_filter(TMX,sgmx,sgmy,npnts,ilm1,ilm2,jlm1,jlm2);


hmsk=HH;
hmsk(HH<0)=nan;
xlim1 = 250;
xlim2 = 1210;
ylim1 = 50;
ylim2 = 1100;


c1=0;
c2=24;
nint=240;
%CMP = colormap_discrete01(c1,c2);
CMP = create_colormap7(nint,c1,c2);
cmp = CMP.colormap;
cint= CMP.intervals;

figure(1); clf;
pcolor(hmsk); shading flat;
colormap([0 0 0]);
freezeColors;

hold on
pcolor(TMXf); shading flat;

colormap(cmp);
caxis([c1 c2]);
axis('equal');
set(gca,'xlim',[xlim1 xlim2],...
	'ylim',[ylim1 ylim2],...
	'Color',[0.5 0.5 0.5]);
set(gca,'xtick',[],'ytick',[]);
clr=[0.9 0.9 0.9];
plot_gridlines(45,10,1,clr,LON,LAT);
stl=sprintf('Timing of Ctr=%3.2f*max(Ctr), %4.1f-%4.1f, GaussF sgm= %4.2f',...
	    rmax,abs(zz1),abs(zz2),sgmx);
title(stl,'Fontsize',12,'Interpreter','none');
ch = colorbar;
set(ch,'ticklength',0.02,...
       'Ticks',[0:2:24],...
       'Fontsize',14);

bottom_text(btx,'pwd',1);

