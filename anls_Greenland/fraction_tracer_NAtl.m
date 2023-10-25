% Plot Tracer fraction in the N.Atl
% monthly-mean, depth-averaged
% tracer concentrations 
% The Data are extracted in extr_trcr_mnth.m
% Option: plot anomalies relative to reference Year/month

%NOTE: the code uses depth-averaged tracer concentration
%and need to be used with caution
%for FW content and dlt S calculation
% layer-averaged tracer cannot be used, instead  
%layer depth-integrated (mass) tracer content has 
% to be used
%
% The Data are extracted in extr_trcr_mnth.m  <---- do not use, saves mean tr
%    output files: trcr_dpthav_*.mat 
% use: extr_MassTrcr_month.m <-- saves depth integrated mass by specified layers
% mat files: MassTr01_lrs_*.mat for Greenland tracer


addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers
startup;

close all
clear

yr1 = 1993;
yr2 = 2016;
nbx = 5; % plot boxes =1,..., nbx

s_fig = 0;
s_mat = 2; % =0 - do not save mat file
           % = 1 save
	   % = 2 skip extraction, only plot output
%dref = datenum(2006,1,15); % reference year/month, = 0 - no reference 
%ilev = 1; %=1 - mixed layer, 0-50; =2 - 50-150m

regn = 'ARCc0.08';
expt = 110;  
pthfig  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/fig_NAtl/',expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
fmout   = sprintf('%sfraction_trcr_NAtl.mat',pthmat);
btxt    = 'fraction_trcr_NAtl.m';

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
f_pltbox=0;
BX = sub_define_boxes(HH,LON,LAT,f_pltbox);
%fbx=sprintf('%sboxes_map',pthfig);
%print('-dpng','-r250',fbx);
% Swap Beaufort sea
%dmm = BX;
%BX(5) = dmm(7);
%BX(6) = dmm(5);
%BX(7) = dmm(6);
%BX = orderfields(BX,[1,2,3,4,7,5,6]);
%keyboard

[XX,YY] = meshgrid((1:nn),(1:mm));
for ib=1:nbx
  iBG = BX(ib).IJ;
  INP = inpolygon(XX,YY,iBG(:,1),iBG(:,2));
  IN = find(INP==1);
  BX(ib).IN = IN;
end


chck=1;
if chck>0
  figure(10); clf;
  contour(HH,[0 0],'k');
  hold on;
  contour(HH,[-5000:1000:-10],'Color',[0.5 0.5 0.5]);
%  contour(HH,[-1000 -1000],'b');
  for ib=1:nbx
    iBG=BX(ib).IJ;
    plot(iBG(:,1),iBG(:,2),'r.-');
  end
%  contour(LON,[-153 -153],'g');
  contour(LON,[-150:30:150],'Color',[0.8 0.8 0.8]);
%  contour(LAT,[70 70],'g');
  contour(LAT,[50:10:80],'Color',[0.9 0.9 0.9]);
  [JJ,II]=ind2sub(size(HH),IN);
  
  for ib=1:nbx
    X=BX(ib).XY(:,1);
    Y=BX(ib).XY(:,2);
    iBG=BX(ib).IJ;
    x0=mean(X);
    y0=mean(Y);
    for kk=1:4
      x=X(kk);
      y=Y(kk);
      dx=x-x0;
      dy=y-y0;
      
      i0=iBG(kk,1)+20*sign(dx);
      j0=iBG(kk,2)+20*sign(dy);
      stl{1}=sprintf('%5.2fN',y);
      stl{2}=sprintf('%5.2fE',x);
      tx=text(i0,j0,stl,'Fontsize',12);
      set(tx,'Color',[0 0 1]);
    end
    
  end
  
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
	for ilv=1:3
	  fprintf('Reading: %i/%2.2i, Tracer %i, Lev %i\n',...
		  iyr,imo,nTr,ilv);
	  Tr = squeeze(TRCR(ilv).TR(nTr,:,:));
	  Tr(Tr<=0)=nan;
  % Integrate Tracer over 1 region of interest per layers
  % To get average Tracer concentration within the region
	  dz = abs(TRCR(ilv).depth_av(2)-TRCR(ilv).depth_av(1));
	  
	  for ib=1:nbx % regions
	    IN=BX(ib).IN;
            volR = sum(Acell(IN)*dz); % total vol, region
      	    mnt = nansum(Tr(IN).*Acell(IN)*dz)/volR; % avrg Tr conc in region
	    BX(ib).mnTR(cc,ilv,nTr) = mnt; % avrg Tr conc in region
	  
  % Integrate Tracer over the whole domain and all depths - total mass
  % rr = mass tracer in the domain / total mass of the tracer
	    intgrTr_dom = nansum(nansum(Tr.*Acell*dz));
	    intgrTr_rgn = nansum(nansum(Tr(IN).*Acell(IN)*dz));
	    rr=intgrTr_rgn/intgrTr_dom; 
	    BX(ib).sumTdom(cc,ilv,nTr) = rr;
	  end
	  
 %if iyr>1995; keyboard; end
	end

      end

    end % month
  end   % year
  if s_mat==1
    fprintf('Saving %s\n',fmout);
    save(fmout,'BX');
  end
else
  fprintf('Loading %s\n',fmout);
  load(fmout);
end
  


%keyboard

% sum over all tracers:
%clear atot
%for ilv=1:2
%  dmm = squeeze(mnTR(:,ilv,:));
%  smm = sum(dmm,2);
%  atot(:,ilv) = smm;
%end

%
% ----------------------
%
% ==========================================================
%   Plot fraction of tracers in the BG wrt to the
% total tracer mass in the domain
% by depth intervals (saved in the extracted time ser.)
% ==========================================================
% Depth level to plot:
ilv = 1; % 0-50m
%ilv = 2; % 50-150
%ilv = 3; % whole depth
switch(ilv)
 case(1)
  smm='0-50m';
 case(2)
  smm='50-150m';
 case(3)
  smm='0-btm';
end

f_pfrc = 1; % plot fraction = Tracer intgr over region/total tracer

if f_pfrc==1
  for ib=1:nbx;
  nm=BX(ib).Name;
  fprintf('Box %i:  %s\n',ib,nm);
  
  sumTdom=BX(ib).sumTdom;
  mnTR=BX(ib).mnTR;

  nrc=length(mnTR);
  yrs=[0:nrc-1]/12+yr1;

  figure(ib); clf;
  % Greeland Runoff
  dmm = squeeze(sumTdom(:,ilv,1));
  axes('position',[0.08 0.71 0.85 0.22]);
  plot(yrs,dmm);
  yl1=0;
  yl2=1.15*max(max(dmm));
  set(gca,'tickdir','out',...
	  'ylim',[yl1 yl2],...
	  'xlim',[yr1 max(yrs)],...
	  'xgrid','on',...
	  'xtick',[yr1:yr2],...
	  'ygrid','on');
  sll=sprintf('%s, Fraction Greenland Trc., %s',nm,smm);
  title(sll);

  % Mackenzie River
  dmm = squeeze(sumTdom(:,ilv,2));
  dmm1 = squeeze(sumTdom(:,ilv,5));
  axes('position',[0.08 0.4 0.85 0.22]);
  plot(yrs,dmm);
  hold
  plot(yrs,dmm1,'r');
  yl1=0;
  yl2=1.15*max(max(dmm1));
  set(gca,'tickdir','out',...
	  'ylim',[yl1 yl2],...
	  'xlim',[yr1 max(yrs)],...
	  'xgrid','on',...
	  'xtick',[yr1:yr2],...
	  'ygrid','on');
  sll=sprintf('%s, Fraction Mkz(b), Bering(r), %s',nm,smm);
  title(sll);

  % Eurasian Rivers
  axes('position',[0.08 0.1 0.85 0.22]);
  dmm = squeeze(sumTdom(:,ilv,3));
  dmm2 = squeeze(sumTdom(:,ilv,4));
  plot(yrs,dmm); % E. Eurasian
  hold
  plot(yrs,dmm2,'r'); % W. Eurasian
  yl1=0;
  yl2=1.15*max(max(dmm2));
  set(gca,'tickdir','out',...
	  'ylim',[yl1 yl2],...
	  'xlim',[yr1 max(yrs)],...
	  'xgrid','on',...
	  'xtick',[yr1:yr2],...
	  'ygrid','on');
  sll=sprintf('%s, Fraction E.(b) W.EurasRiv, %s',nm,smm);
  title(sll);

  bottom_text(btxt,'pwd',1);

  if s_fig==1
    fgnm=sprintf('%sfrTrc_%s_Lev%i',pthfig,nm,ilv);
    fprintf('Saving %s\n',fgnm);
    print('-dpng','-r250',fgnm);
  end

  end
end

% ==============================
%
% Fraction in different regions
% for different depths:
%
% ==============================
for ilv=1:3
  for ib=1:nbx;
    nm=BX(ib).Name;
    if ilv==1, xlb{ib}=nm; end;
    fprintf('Box %s\n',nm);
  
    sumTdom=BX(ib).sumTdom;
    mnTR=BX(ib).mnTR;

  % Greeland Runoff
    dmm = squeeze(sumTdom(:,ilv,1));
    
    A(ilv,ib)=mean(dmm(end-11:end));
  end
end
figure(30); clf;
axes('Position',[0.08 0.45 0.88 0.5]);
bar(A',0.95);
sll=sprintf('Fraction of Greenland Tracer in 2016');
title(sll);
set(gca,'xticklabel',xlb,...
	'tickdir','out');
lg=legend('0-50','50-150','btm');

btx='fraction_tracer_NAtl.m';
bottom_text(btx,'pwd',1,'Position',[0.08 0.3 0.5 0.1])

% ==================
% Mean concentration
% of individual tracers 
% inside the boxes
% =================
for ib=1:nbx;
  nm=BX(ib).Name;
  mnTR=BX(ib).mnTR;

  fprintf('Box %s\n',nm);

  nrc=length(mnTR);
  yrs=[0:nrc-1]/12+yr1;

  figure(10+ib); clf;
%  axes('position',[0.08 0.7 0.85 0.25]);
  axes('position',[0.08 0.71 0.85 0.22]);
  dmm = squeeze(mnTR(:,ilv,1)); %
  plot(yrs,dmm,'linewidth',2);
  hold on
%  plot(yrs,dmm5,'k','linewidth',2); % Bering Strait water

  yl1=0.95*(min(dmm));
  yl2=1.1*max(dmm);
  set(gca,'tickdir','out',...
	  'ylim',[yl1 yl2],...
	  'xlim',[yr1 max(yrs)],...
	  'xtick',[yr1:yr2],...
	  'xgrid','on',...
	  'ygrid','on',...
	  'xminortick','on',...
	  'yminortick','on');
  sll=sprintf('%s, Conc GreenTr, %s',nm,smm);
  title(sll);

% Pacific Water/Bering str.   
  axes('position',[0.08 0.4 0.85 0.22]);
  dmm5 = squeeze(mnTR(:,ilv,5)); %
  plot(yrs,dmm5,'linewidth',2); % Bering Strait water
  yl1=0.95*(min(dmm5));
  yl2=1.1*max(dmm5);
  set(gca,'tickdir','out',...
	  'ylim',[yl1 yl2],...
	  'xlim',[yr1 max(yrs)],...
	  'xtick',[yr1:yr2],...
	  'xgrid','on',...
	  'ygrid','on',...
	  'xminortick','on',...
	  'yminortick','on');
  sll=sprintf('%s, Conc Bering, %s',nm,smm);
  title(sll);
  
  
  
  dmm2 = squeeze(mnTR(:,ilv,2)); %
  dmm3 = squeeze(mnTR(:,ilv,3)); %
  dmm4 = squeeze(mnTR(:,ilv,4)); %
%  axes('position',[0.08 0.35 0.85 0.25]);
  axes('position',[0.08 0.1 0.85 0.22]);
  dmm = squeeze(mnTR(:,ilv,4)); % W. Eurasian Rivers
  plot(yrs,dmm2,'linewidth',2);
  hold on
  plot(yrs,dmm3,'r','linewidth',2);
  plot(yrs,dmm4,'g','linewidth',2);
  %plot(yrs,dmm5,'k','linewidth',2);

  yl1=0.95*(min(dmm3));
  yl2=1.1*max([max(dmm3), max(dmm4),max(dmm2)]);
  set(gca,'tickdir','out',...
	  'ylim',[yl1 yl2],...
	  'xlim',[yr1 max(yrs)],...
	  'xtick',[yr1:yr2],...
	  'xgrid','on',...
	  'ygrid','on',...
	  'xminortick','on',...
	  'yminortick','on');
  sll=sprintf('%s, Conc MkzTr(b), E.(r) & W.(g)Riv. Tracers, %s',nm,smm);
  title(sll);

  btxt='fraction_trcr_NAtl.m';
  bottom_text(btxt,'pwd',1);


%keyboard
  if s_fig==1
    fgnm=sprintf('%sconcTrcrs_%s_Lev%i',pthfig,nm,ilv);
    fprintf('Saving %s\n',fgnm);
    print('-dpng','-r250',fgnm);
  end

end

