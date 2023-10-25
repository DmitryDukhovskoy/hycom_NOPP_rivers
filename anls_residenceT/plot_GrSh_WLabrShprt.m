% Combine all experiments with particles:
% N experiments with particles
% released at different depth levels
%  
% Plot particles released at SW Greenland shelf and 
% NW Labr shelf 
%
% First, need to combine all experiments/runs 
% For Gr Shelf: fget = 1 run plot_GrShprt_paths.m
% for NWLabr SHelf:      plot_WLabrShprt_paths.m
%  
% Compare different levels
%
% 
addpath /Net/Movies0/ddmitry/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /home/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers;
%addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
%addpath /Net/Movies0/ddmitry/MyMatlab/seawater
startup

close all
clear

fget=0;  % =0 - load previously loaded data

cc = 0;
ilr = 3;  % layer to plot

NSM_Gr=[1,2,3,4,5,6]; % each epxeriment released 100 particles at 1 depth
NSM_WL=[1,2,3,4]; % WLabr only 3 runs for each layer

%NSM=[1,2,3,4];
ngrm=length(NSM_Gr);
nwlm=length(NSM_WL);

SLR=[10,15,23,31]; % depth levels of the particels
ZLR=[50,90,150,450];  % nominal depths of particles in deep ocean
nlr=length(SLR);

%YRPLT=[1994 2002 2019];
YRPLT = 1994;
nplt = length(YRPLT);


%slr = 15; % vertical layer
%s_par=1; % parallel session 
%f_plt=1; % plot prtcles

s_fig = 1;

fprintf('Save fig %i, Layer=%i\n',s_fig,ilr);

%cc = 71; % last SAVED frame - if need to restart from frame # 
      % the next frame will be cc+1
      % = 0 - start from beginning

rg = 9806; 

regn = 'ARCc0.08';
expt = 110;  
%pthfig  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/fig_trac/',expt);
pthfig  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/%s/%3.3i/fig_GrGprtcl/',regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/110/data_GG_prt/';
pthmat2 = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/110/data_WLabr_prt/';

txtb = 'plot_GrSh_WLabrShprt.m';

% ------------------------
% TOPO
% ------------------------
ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

hmsk=HH;
hmsk(HH<0)=nan;

xlim1 = 220;
xlim2 = 1250;
ylim1 = 30;
ylim2 = 1100;

% layer interface depths:
% in N Atl:
zz=[                     0
                        -1
         -2.80000007967061
         -6.04000002390118
         -10.7199998326917
         -15.6499996414823
         -21.4599995777458
         -28.3299994502728
         -36.3299994502728
         -44.3299994502728
         -52.3299994502728
         -60.3299994502728
         -68.3299994502728
         -76.3299994502728
         -84.3299994502728
         -92.3299994502728
         -100.329999450273
         -108.329999450273
         -116.329999450273
         -124.329999450273
         -132.329999450273
         -140.329999450273
         -148.329999450273
         -156.329999450273
         -166.329999450273
         -182.730000087638
         -218.650001234894
         -261.030001362367
         -311.050001872259
         -370.070002382151
         -439.709999577746
         -521.889997793124
         -642.060033995449
         -833.055521452108
         -1631.30343089531
         -2130.04884186818
          -2911.1784563899
         -3125.48785879659
         -3385.88090896996
         -3385.88090896996
         -3385.88090896996
         -3385.88090896996];

f_dps = 1;
if f_dps==1
% Read layer thicknesses:
  pthbin = '/nexsan/hycom/ARCc0.08_112/data/1993/';
  fina = sprintf('%s112_archm.1993_002_12.a',pthbin);
  finb = sprintf('%s112_archm.1993_002_12.b',pthbin);

  [ZM,ZZ] = sub_zz_zm(fina, finb,HH);
% on the shelf depth of the layer:
  ZM(31,539,529)
  
  jz=10;
  jz=15;
  jz=23;
  jz=31;
  d1=squeeze(ZZ(jz-1,:,:));
  d2=squeeze(ZZ(jz+1,:,:));
  I=find(HH<-1000);
  md1=nanmean(d1(I));
  md2=nanmean(d2(I));
  fprintf('Layer %i, top/btm zz=%8.4f/%8.4f \n',jz,md1,md2);

  ii1=518;
  ii2=557;
  jj1=470;
  jj2=470;
 
  zzx=squeeze(ZZ(:,jj1,ii1:ii2));
  hb=HH(jj1,ii1:ii2);

  figure(20); clf;
  hold on;
  plot(hb,'k-','Linewidth',2);  

  [a1,a2]=size(zzx);
  for ii=1:a1
    zl=zzx(ii,:);
    plot(zl,'-');
  end; 
  set(gca,'xlim',[15 40],...
          'ylim',[-3500 0]);

keyboard
end


fnmout=sprintf('%sGrShprt_comb.mat',pthmat);
fnmout2=sprintf('%sWLabrShprt_comb.mat',pthmat2);

% Combine all experiments with Lagrangian particles together 
% for missing dates - use the previous number and locations of particles
PPG = sub_GrSh(fnmout);
PPL = sub_GrSh(fnmout2);


Hmsk=HH*0;
Hmsk(HH<0)=1;

% Region of interest - double check
% with dS_FWC_timeseries_SubpolarGyre.m
fspg='SPG_noNorth_indx.mat';  % get IGR indices of the region
load(fspg);
IGR(end+1,:)=IGR(1,:);



% Particles release location:
%fmatLP = sprintf('%sGrSh_part-lr10_1993-01.mat',pthmatLP);
%Ip = PRTCL.TRACK(1).I;
%Jp = PRTCL.TRACK(1).J;
slr=[];
fina=[];
finb=[];
Np0=1e4;
nsim=1;
[bmm,dmm,amm]=sub_seed_GrSh_prt(slr,HH,LON,LAT,fina,finb,Np0,nsim);
Ip = bmm.TRACK.I;
Jp = bmm.TRACK.J;

%
% Labrador Sea:
[bmm,dmm,amm]=sub_seed_WLabr_prt(slr,HH,LON,LAT,fina,finb,Np0,nsim);
Ilb = bmm.TRACK.I;
Jlb = bmm.TRACK.J;


%cc=0;
figure(1);
set(gcf,'Position',[534 76 1195 1266]);
%set(gcf,'Position',[9 5 707 797]);

%for itt=733:733
for iyr=1:nplt
  figure(1); clf;
  pcolor(Hmsk); shading flat;
  colormap([0 0 0; 1 1 1]);
  freezeColors;
  caxis([0 1]);

  hold on;
  
  contour(HH,[-8000:1000:-100],'Linewidth',1,'Color',[0.8 0.8 0.8]);
  contour(HH,[-1000 -1000],'Linewidth',1.6,'Color',[0.5 0.5 0.5]);
  contour(HH,[-500 -500],'Linewidth',1,'Color',[0.8 0.8 0.8]);

% Plot Lagr Prt WLabr Shelf Particles
% release locations
plot(Ip,Jp,'.','Color',[0.7 0.8 1]);
plot(Ilb,Jlb,'.','Color',[1 0.8 0.7]);

% SPNA domain polygon:
  plot(IGR(:,1),IGR(:,2),'-','Linewidth',1.6,...
       'Color',[0.6 0. 0.]);

% Plot GrSh particles:
  TM  = PPG(1).TM;
  DV   = datevec(TM);
  nrc  = length(TM);

  yr0=YRPLT(iyr);
  I=find(DV(:,1)== yr0);
  itt=I(end);
  dv1=datevec(TM(itt));
  fprintf('Plotting %i/%i/%i\n',dv1(1:3));
  
  
%keyboard
		lr=SLR(ilr);

  XP = PPG(1).XPcmb;
  YP = PPG(1).YPcmb;
  ZL = PPG(1).ZLcmb;
		IZ=find(ZL==lr);
  clrg=[0 0.4 0.8];

		plot(XP(IZ,itt),YP(IZ,itt),'.',...
					'Markersize',31,...
					'Color',clrg);

% Plot WLabrSh particles:
  TM  = PPL(1).TM;
  DV   = datevec(TM);
  nrc  = length(TM);

  I=find(DV(:,1)== yr0);
  itt=I(end);
  dv1=datevec(TM(itt));

  fprintf('WLabr Sh: Plotting %i/%i/%i\n',dv1(1:3));

  XPl = PPL(1).XPcmb;
  YPl = PPL(1).YPcmb;
  ZLl = PPL(1).ZLcmb;
		IZl=find(ZLl==lr);

  clrl = [0.8 0.3 0];
		plot(XPl(IZl,itt),YPl(IZl,itt),'.',...
					'Markersize',26,...
					'Color',clrl);

  axis('equal');
  set(gca,'xlim',[xlim1 xlim2],...
	  'ylim',[ylim1 ylim2]);
  set(gca,'xtick',[],'ytick',[]);
    
  stt = sprintf('%2.2i/%2.2i/%i, Lr=%i, Z=%i m',dv1(3),dv1(2),dv1(1),SLR(ilr),ZLR(ilr));
  title(stt);
 
% Legend
  axes('Position',[0.15 0.06 0.2 0.05]);
  hold on;
  plot(0.5,1,'.','Markersize',28,'Color',clrg);
  text(0.7,1,'West Gr Sh','Fontsize',14);
  plot(0.5,0.5,'.','Markersize',28,'Color',clrl);
  text(0.7,0.5,'North Labr Sh','Fontsize',14);
  set(gca,'xlim',[0.45 1.3],...
          'ylim',[0.4 1.2],...
          'visible','off')



  bottom_text(txtb,'pwd',1,'Fontsize',7);
%    drawnow
  cc = cc+1;
%keyboard    
  if s_fig>0
    fgnm = sprintf('%s008-110_GrWLabrSh_prtcl_lr%2.2i-%4.4i%2.2i%2.2i',pthfig,ilr,dv1(1:3));

    fprintf('Saving figure %s\n',fgnm);
    print('-dpng','-r200',fgnm);
  end
  %fprintf('1 record: %8.5f min\n\n',toc/60);
		     %    keyboard

end







