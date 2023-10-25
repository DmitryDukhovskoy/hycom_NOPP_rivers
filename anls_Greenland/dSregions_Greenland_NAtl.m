% Estimate S change in the upper layer due to surplus FW flux
% from Greenland 
% Use estimate of Greenland runoff fraction
% for specified boxes
% the Fraction is calculated in fraction_tracer_NAtl.m
%NOTE: the code uses depth-averaged tracer concentration
% extracted in extr_tracer_month.m -> fraction_tracer_NAtl.m
%and need to be used with caution
%for FW content and dlt S calculation
% layer-averaged tracer cannot be used, instead  
%layer depth-integrated (mass) tracer cnotent has 
% to be used
% The Data are extracted in extr_trcr_mnth.m  <---- do not use, saves mean tr
% use: extr_MassTrcr_month.m <-- saves depth integrated mass by specified layers
% mat files: MassTr01_lrs_*.mat for Greenland tracer
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers
startup;

close all
clear

s_fig = 0;

nbx = 5; % # of regions
regn = 'ARCc0.08';
expt = 110;  
pthfig  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/fig_NAtl/',expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
fmout   = sprintf('%sfraction_trcr_NAtl.mat',pthmat);

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

% Define Regions:
BX = sub_define_boxes(HH,LON,LAT,0);
%fbx=sprintf('%sboxes_map',pthfig);
%print('-dpng','-r250',fbx);



[XX,YY] = meshgrid((1:nn),(1:mm));
for ib=1:nbx
  iBG = BX(ib).IJ;
  INP = inpolygon(XX,YY,iBG(:,1),iBG(:,2));
  IN = find(INP==1);
  BX(ib).IN = IN;
end


chck=0;
if chck>0
  figure(10); clf;
  contour(HH,[0 0],'k');
  hold on;
  contour(HH,[-5000:1000:-10],'b');
  for ib=1:nbx
    iBG=BX(ib).IJ;
    plot(iBG(:,1),iBG(:,2),'r.-');
  end
  contour(LON,[-153 -153],'g');
  contour(LON,[-130 -130],'g');
  contour(LAT,[73 73],'g');
  contour(LAT,[83 83],'g');
  [JJ,II]=ind2sub(size(HH),IN);
  plot(II,JJ,'y.');
  axis('equal');
end


fprintf('Loading %s\n',fmout);
load(fmout);


% sum over all tracers:
%clear atot
%for ilv=1:2
%  dmm = squeeze(mnTR(:,ilv,:));
%  smm = sum(dmm,2);
%  atot(:,ilv) = smm;
%end

[dGR,YGR]=sub_read_Greenland; % Greenland anom, km3/yr


% Get S fields from HYCOM
yr = 1993;
mo = 1;
dnmb = datenum(yr,mo,15);
iday = dnmb-datenum(yr,1,1)+1;
fprintf('Extracting %i/%i\n',yr,mo);
expt=110;

pthbin  = sprintf('/nexsan/archive/ARCc0.08_110/data/%4.4i/',yr);
fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);
fld='salin';
[SS,n,m,l] = read_hycom(fina,finb,fld,'r_layer',2);
SS=squeeze(SS);
SS(SS>1e20)=nan;


% ==========================================================
%   Derive fraction of tracers in the region wrt to the
% total tracer mass in the domain
% by depth intervals (saved in the extracted time ser.)
% Only N.Atl tracers - estimate only Greenland FW
% ==========================================================
yr1=1993; % start of the simulation
ilv=1;  % Mixed Layer, 0-50 m
dz=50;  %m,  need to change if use different layers
itr=1;  % Greenland Tracer
Frtr=squeeze(BX(1).sumTdom(:,ilv,itr)); % fraction of tracer in the region

nrc=length(Frtr);
yrs=[0:nrc-1]/12+yr1;

% Reconstruct surplus Greenland FW flux, km3
iyr = find(YGR==yr1);
clear gflx
gflx(1,1) = dGR(1)/12;
clear it ib
for it=2:nrc
  yr0=floor(yrs(it));
  mo=round((yrs(it)-yr0)*12)+1;
  iyr = find(YGR==yr0);
  gflx(it)=gflx(it-1)+dGR(iyr)/12; % km3 of surplus FW
end
gflx=gflx(:);

Need to use Bamber 2018 river runoff:
Modify code as follows
frv = sprintf('%sGreenland_cumFWFlux.mat',pthmat);
fprintf('f_griv %i, Loading %s\n',f_griv,frv);
load(frv);
gflx = cFWF; % km3 

f_pgr=0; % plot cumulative anomaly Greenland FW flux
if f_pgr==1
  figure(5); clf;
  axes('Position',[0.08 0.5 0.8 0.43]);
  plot(yrs,gflx,'Linewidth',2); % 
  set(gca,'tickdir','out',...
	  'ylim',[0 1.1*max(gflx)],...
	  'xlim',[1993 2016],...
	  'xgrid','on',...
	  'xtick',[yr1:yr2],...
	  'xminortick','on',...
	  'yminortick','on',...
	  'ygrid','on');
  sll=sprintf('Cumulative Surplus Greenland FWF, km3');
  title(sll);

  
  btxt='dSregions_Greenland_NAtl.m';
  bottom_text(btxt,'pwd',1,'position',[0.02 0.3 0.9 0.1]);

  if s_fig==1
    fgnm=sprintf('%sCumulativeGreenFWF',pthfig);
    fprintf('Saving %s\n',fgnm);
    print('-dpng','-r250',fgnm);
  end
  
end


% Need to extract monthly mean S
% for layers

clear GFWbx* dS*
for ib=1:nbx;
  nm=BX(ib).Name;
  Frtr=squeeze(BX(ib).sumTdom(:,ilv,itr)); % fraction of tracer in the region
  
% Volume of surplus Greenland FW in the boxes, km3
% convert to m of FW, i.e. normalize by area of region
  IN = BX(ib).IN;
  Arg=sum(Acell(IN)); % m2
  Sav=nanmean(SS(IN)); % mean S
  Vol = Arg*dz;  % region volume, m3
  GFWbx(:,ib)=gflx.*Frtr.*1e9; % m3   
  GFWbxM(:,ib)=gflx.*Frtr.*1e9/Arg; % m3 -> m
% Estimate salinity changed dltS = S-dS:
  dS(:,ib)=Vol*Sav./(Vol+GFWbx(:,ib));
  
end


% =================================================
% Plot FW content change due to surplus Greenland FW
% Fluxes, m
%
% =================================================
POS=[0.08 0.75 0.88 0.15;...
     0.08 0.53 0.88 0.15;...
     0.08 0.31 0.88 0.15;...
     0.08 0.08 0.88 0.15];
     
yr2=2016;

figure(1); clf;
for ibb=1:nbx;
  ib=ibb;
  if ibb>4
    figure(11); clf;
    ib=ibb-4;
  end
  
  nm=BX(ibb).Name;
  gfwf=GFWbxM(:,ibb);
  
  nrc=length(gfwf);
  yrs=[0:nrc-1]/12+yr1;

  pos = POS(ib,:);
  axes('position',pos);
  plot(yrs,gfwf,'Linewidth',2);
  yl1=0;
  yl2=1.15*max(max(gfwf));
  set(gca,'tickdir','out',...
	  'ylim',[yl1 yl2],...
	  'xlim',[yr1 max(yrs)],...
	  'xgrid','on',...
	  'xtick',[yr1:yr2],...
	  'yminortick','on',...
	  'ygrid','on',...
	  'fontsize',12);
  sll=sprintf('Surplus Greenl.FW (m) in %s, 0-50m',nm);
  title(sll);

end

btxt='dSregions_Greenland_NAtl.m';
bottom_text(btxt,'pwd',1);

if s_fig==1
  fgnm=sprintf('%sGreenFWanom_%s_Lev%i',pthfig,nm,ilv);
  fprintf('Saving %s\n',fgnm);
  print('-dpng','-r250',fgnm);
end

figure(2); clf;
for ibb=1:nbx;
  ib=ibb;
  if ibb>4
    figure(12); clf;
    ib=ibb-4;
  end
  
  nm=BX(ibb).Name;
  dmm=dS(:,ibb);
  
  nrc=length(gfwf);
  yrs=[0:nrc-1]/12+yr1;

  pos = POS(ib,:);
  axes('position',pos);
  plot(yrs,dmm,'Linewidth',2);
  yl2=1.0001*max(max(dmm));
  yl1=0.9999*min(min(dmm));
  set(gca,'tickdir','out',...
	  'ylim',[yl1 yl2],...
	  'xlim',[yr1 max(yrs)],...
	  'xgrid','on',...
	  'xtick',[yr1:yr2],...
	  'xminortick','on',...
	  'yminortick','on',...
	  'ygrid','on');
  sll=sprintf('dlt(S) due to Surplus Greenl.FW,  %s, 0-50m',nm);
  title(sll);

end

btxt='dSregions_Greenland_NAtl.m';
bottom_text(btxt,'pwd',1);

if s_fig==1
  fgnm=sprintf('%sSchange_GreenFWanom_%s_Lev%i',pthfig,nm,ilv);
  fprintf('Saving %s\n',fgnm);
  print('-dpng','-r250',fgnm);
end




