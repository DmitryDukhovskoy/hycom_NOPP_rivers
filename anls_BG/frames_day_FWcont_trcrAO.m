% Plot monthly-mean, depth-integrated 
% similar to
% daily tracer conc fraction but scaled
% to represent fwc (cm) within the layer
% The Data are extracted in extr_trcr_day.m
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers
startup;

close all
clear

s_fig = 1;
ilv = 1;  % level to plot/analyze
cc=0; % last SAVED frame - if need to restart from frame # 
      % the next frame will be cc+1
      % = 0 - start from beginning

% If need to continue from the last saved frame NN
% fr_last = NN, 
% the code will go from year yr1 and finds 
% what year/day it needs to continue
% make sure that yr1 correspond to frame 0001 
fr_last = 0; % last saved frame - if need to restart from frame # 
        % = 0 - start from beginning
yr1 = 1993;
yr2 = 2016;

fprintf('Last saved frame %i in %i-%i\n',fr_last,yr1,yr2);


regn = 'ARCc0.08';
expt = 110;  
pthfig=sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/frames_day_trcrAO/',...
		  expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
txtb='frames_day_FWcont_trcrAO.m';

ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);
ACell=DX.*DY;

xlim1 = 20;
xlim2 = nn;
ylim1 = 500;
ylim2 = 2000;

hmsk=HH;
hmsk(HH<0)=nan;

% Load river runoff 
% to scale Tr Conc to FW content
f_getrv=0;
  frv=sprintf('%strcr_rivers_fwflux.mat',pthmat);

if f_getrv==1
  mday = [31;28;31;30;31;30;31;31;30;31;30;31];
  BB.hycom_topo   = HH;
  BB.hycom_LON    = LON;
  BB.hyom_LAT     = LAT;
  BB.area_cf      = 16e6;  % 16e6 - mean cell area, m2, for runoff m/s->m3/sec
  BB.Area_cell_m2 = ACell; % if empty, - use area_cf;
  BB.river_flag   = 'amer';
  RR = river_FWflux_byTrcrs(BB);
  FWF(2).Name='Mackz';
  FWF(2).FWF_km3_mo=RR.FWF_m3_mo*1e-9; % km3/mo
  FWF(2).TM=RR.TM;
  FWF(2).cumFW_km3=cumsum(RR.FWF_m3_mo*1e-9); 
  
  BB.river_flag   = 'W-euras';
  RR = river_FWflux_byTrcrs(BB);
  FWF(3).Name='W.Eurs';
  FWF(3).FWF_km3_mo=RR.FWF_m3_mo*1e-9; % km3/mo
  FWF(3).TM=RR.TM;
  FWF(3).cumFW_km3=cumsum(RR.FWF_m3_mo*1e-9); 
  
  BB.river_flag   = 'E-euras';
  RR = river_FWflux_byTrcrs(BB);
  FWF(4).Name='E.Eurs';
  FWF(4).FWF_km3_mo=RR.FWF_m3_mo*1e-9; % km3/mo
  FWF(4).TM=RR.TM;
  FWF(4).cumFW_km3=cumsum(RR.FWF_m3_mo*1e-9); 

  cc=0;
  for iyr=1993:2016
    for im=1:12
      cc=cc+1;
      bfwf(cc)=2779/365*mday(im); % use HYCOM FWflux Bering
    end
  end
  FWF(5).Name='Bering';
  FWF(5).FWF_km3_mo=bfwf; % km3/mo
  FWF(5).TM=RR.TM;
  FWF(5).cumFW_km3=cumsum(bfwf); 

  save(frv,'FWF');
  
  
else
  fprintf('Loading rivers\n',frv);
  load(frv);
  
end

fplt_cumR=0;
if fplt_cumR==1
  YRS=[1993:1/12:2016.99];
  figure(11); clf;
  axes('Position',[0.08 0.5 0.85 0.4]);
  hold on
  clear ltx
  for ik=2:5
    plot(YRS,FWF(ik).cumFW_km3);
    fmn=mean(FWF(ik).FWF_km3_mo)*12; %km3/yr
    ltx{ik-1}=sprintf('%s: %6.1f km3',...
		    FWF(ik).Name,fmn);
  end
  lpp=legend('Mackz','W.Eurs','E.Eurs','Bering');
  set(lpp,'Position',[0.11 0.74 0.09 0.12]);
  text(1993.5, 25000,ltx,'Fontsize',11);
  
  
  set(gca,'tickdir','out',...
	  'xlim',[1993 2016],...
	  'ylim',[0 7.2e4],...
	  'xtick',[1990:2:2020],...
	  'ytick',[0:10000:80000],...
	  'xgrid','on',...
	  'ygrid','on',...
	  'Fontsize',14);
  title('Cumul. FWFlx (km3), Tracer Tracekd Rivers and BeringStr.');
  
  bottom_text(txtb,'pwd',1,'position',[0.04 0.4 0.5 0.02]);
  
  
end



nint = 320;
c1 = 0; %
c2 = 100; %
%CMP = create_colormap2_3(nint,c1,c2);
CMP = colormap_sclr2(nint,c1,c2);
cmp = CMP.colormap;
cmp(1,:) = [1 1 1];
cmp(2,:) = [1 1 1];
cmp(3,:) = [1 1 1];
cmp(4,:) = [1 1 1];
cmp(5,:) = [1 1 1];
cmp(6,:) = [1 1 1];
cmp(7,:) = [1 1 1];
cmp(8,:) = [1 1 1];
cmp(9,:) = [1 1 1];
cmp(10,:) = [1 1 1];
cmp = smooth_colormap(cmp,18);
cmp = smooth_colormap(cmp,18);
cmp(1,:) = [1 1 1];
cnt = CMP.intervals;

ck = 0;
YRPLT = [yr1:yr2];
nrc=length(YRPLT);

szd=0.44;
szu=0.475;
POS(2,:) = [0.02 0.51 szd szu];
POS(3,:) = [0.465 0.51 szd szu];
POS(4,:) = [0.02 0.06 szd szu];
POS(5,:) = [0.465 0.06 szd szu];


close all
ff=figure('visible','off');


ic1 = 1;
cc=0; 
for ic=ic1:nrc
  iyr = YRPLT(ic);
  fmat = sprintf('%strcr_dpthav_daily_%4.4i.mat',pthmat,iyr);
  
  if ~exist(fmat,'file')
    fprintf('File is missing %s\n',fmat);
    continue
  end

  fprintf('Loading %s\n',fmat);
  load(fmat);
  

  TM = TRCR(1).TM;
  ndy = length(TM);
  for idd = 1:ndy
    tic;
    dnmb = TM(idd);
    dv = datevec(dnmb);
    iyr=dv(1);
    imo=dv(2);
    cc = cc+1;

    if cc<=fr_last,
      fprintf('Frame %4.4i need %4.4i, skipping ...\n',cc,fr_last);
      continue;
    end
    
    clf;
    for nTr = 2:5
      TMR=FWF(nTr).TM;
      itm = find(TMR==datenum(iyr,imo,15));
      cFWF= FWF(nTr).cumFW_km3(itm)*1e9; % m3

      fprintf('Plotting: %2.2i/%2.2i/%2.2i, Tracer %i, Level=%i\n',...
	    dv(1:3),nTr,ilv);
      
      if ilv==1
	dz=abs(TRCR(nTr).depth_av1(2)-TRCR(nTr).depth_av1(1));
        Tr = double(squeeze(TRCR(nTr).TR_lr1(idd,:,:)));
	rmm = double(squeeze(TRCR(nTr).TR_lr2(idd,:,:)));
	dzm=abs(TRCR(nTr).depth_av2(2)-TRCR(nTr).depth_av2(1));
      else
	dz=abs(TRCR(nTr).depth_av2(2)-TRCR(nTr).depth_av2(1));
        Tr = double(squeeze(TRCR(nTr).TR_lr2(idd,:,:)));
	rmm = double(squeeze(TRCR(nTr).TR_lr1(idd,:,:)));
	dzm=abs(TRCR(nTr).depth_av1(2)-TRCR(nTr).depth_av1(1));
      end
      
% Combin all tracer - mass 
% Note that whole depth is missing in 
% the extracted files,
% use only 0-150 m mass
% Will need to change this
      Tr(Tr<=0)=0;
      rmm(rmm<=0)=0;
      Mtr=Tr.*ACell*dz;
      Mtr2=rmm.*ACell*dzm;
      Mall=Mtr+Mtr2; 
%	lTr=log(Tr);
%	lTr(isnan(lTr))=-1000;
%	lTr(HH>=0)=nan;
      Mall(Mall==0)=nan;
      rr=Mtr./nansum(nansum(Mall));
      rr(rr<=0)=nan;
% Scale rr to represent the FWC	within 10 m of water
      lrr = rr*cFWF./ACell*100; % cm of FWF within the layer
      lrr = lrr./dz*10; 
%        lrr=log(rr);
      lrr(isnan(lrr))=-1e30;
      lrr(HH>=0)=nan;


      if ilv==1
	lvl='0-50m';
      else
	lvl='50-150m';
      end

      pst = POS(nTr,:);
      axes('Position',pst);
%	pcolor(lTr); shading flat;
      pcolor(lrr); shading flat;
      hold on;
      contour(HH,[0 0],'k','Linewidth',1);
      contour(HH,[-1000 -1000],'Color',[0.7 0.7 0.7],'Linewidth',0.8);
      caxis([c1 c2]);
      colormap(cmp);
      axis('equal');
      set(gca,'xlim',[xlim1 xlim2],...
	      'ylim',[ylim1 ylim2],...
	      'Color',[0. 0. 0.]);
      set(gca,'xtick',[],'ytick',[]);
      clr=[0.9 0.9 0.9];
      plot_gridlines(90,20,0.5,clr,LON,LAT);

      clear stl
      stl{1} = sprintf('ARCc0.08-%3.3i',expt);
      stl{2} = sprintf('%2.2i/%2.2i/%2.2i',dv(1:3));
      stl{3} = sprintf('Trcr# %i, %s',nTr,lvl);
      tp=text(50,1750,stl);
      set(tp,'Fontsize',12,'Color',[0 0.8 1]);
      title('FWC (cm in 10m)');
%	title(stl,'Fontsize',8);

      if nTr == 3 
	clb=colorbar;
	set(clb,'TickLength',0.02,...
	      'Position',[0.92, 0.55, 0.02, 0.4],...
	      'Fontsize',12);
      end	  
      if nTr == 5 
	clb=colorbar;
	set(clb,'TickLength',0.02,...
	      'Position',[0.92, 0.08, 0.02, 0.4],...
	      'Fontsize',12);
      end	  

    end  % nTr
    set(gcf,'Position',[1199 81 1195 1250]);
% keyboard   
    if s_fig>0
      bottom_text(txtb,'pwd',1,...
		  'position',[0.04 0.05 0.4 0.03],...
		  'Fontsize',8);
      fgnm = sprintf('%sday_FWtrcrLev%2.2i_%5.5i',pthfig,ilv,cc);
      fprintf('Saving figure %s\n',fgnm);
%      print('-djpeg','-r300',fgnm);
      print('-dpng','-r200',fgnm);
    end
    fprintf('1 record: %8.5f min\n\n',toc/60);

  end; % days
end   % year loop





