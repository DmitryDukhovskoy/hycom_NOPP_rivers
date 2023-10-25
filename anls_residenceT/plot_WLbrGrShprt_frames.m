% Plot WLabr and Gr Shelf together
% Combine all experiments with particles:
% N experiments with particles
% released at different depth levels
% 
% Particles are not added during the simulation
% All N prticles seeded at once at initial state
% For each particle: 
% no T or S is tracked, only time and location
% 
addpath /Net/Movies0/ddmitry/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /home/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers;
%addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
%addpath /Net/Movies0/ddmitry/MyMatlab/seawater
startup

close all
clear

shelf1='GrSh';
shelf2='WLabr';

s_fig = 2019; % s_fig = 2003; start from last saved in 2003, frame counter is from it1
iStrt = 2019; 
iEnd  = 2019;
dtt   = 4;    % skip dtt output fields (ndays = dtt*dt, dt = 2 days)


cc = 0;
%dtpth = 90; % plot trajectory over the last n days


NSM=[1,2,3,4,5,6]; % each epxeriment released 100 particles at 1 depth
%NSM=1;
nnm=length(NSM);

%SLR=10;
SLR=[10,15,23,31]; % depth levels of the particels
ZLR=[50,90,150,450];  % nominal depths of particles in deep ocean
nlr=length(SLR);


fprintf('Plot frames with particles: Save fig %i \n',s_fig);

%cc = 71; % last SAVED frame - if need to restart from frame # 
      % the next frame will be cc+1
      % = 0 - start from beginning

rg = 9806; 

regn = 'ARCc0.08';
expt = 110;  
%pthfig  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/fig_trac/',expt);
pthfig = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/%s/%3.3i/frames_GrGprtcl/',regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat1  = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/110/data_GG_prt/';
pthmat2  = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/110/data_WLabr_prt/';
txtb = 'plot_WLbrGrShprt_frames.m';

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

f_dps = 0;
if f_dps==1
% Read layer thicknesses:
  pthbin = '/nexsan/hycom/ARCc0.08_112/data/1993/';
  fina = sprintf('%s112_archm.1993_002_12.a',pthbin);
  finb = sprintf('%s112_archm.1993_002_12.b',pthbin);

  [ZM,ZZ] = sub_zz_zm(fina, finb,HH);
% on the shelf depth of the layer:
 ZM(31,539,529)
 
end


fnmout1=sprintf('%sGrShprt_comb.mat',pthmat1);
fprintf('Loading %s\n',fnmout1);
load(fnmout1);
PPG = PP;

fnmout2=sprintf('%sWLabrShprt_comb.mat',pthmat2);
fprintf('Loading %s\n',fnmout2);
load(fnmout2);
PPL = PP;

PRTG = sub_combPrt(PPG);
PRTL = sub_combPrt(PPL);



nn1=600;
cc1=1;
cc2=nn1;
CMP = create_colormap_paleB(nn1,cc1,cc2);
%CMP = create_colormap6(nn1,cc1,cc2);
CLRF = CMP.colormap;


Hmsk=HH*0;
Hmsk(HH<0)=1;

% Region of interest - double check
% with dS_FWC_timeseries_SubpolarGyre.m
fspg='SPG_noNorth_indx.mat';  % get IGR indices of the region
load(fspg);
IGR(end+1,:)=IGR(1,:);


[XX,YY] = meshgrid((1:nn),(1:mm));

INP = inpolygon(XX,YY,IGR(:,1),IGR(:,2));
INGr= find(INP==1 & HH<0);
II = find(HH<0);


dx=0.45;
POS = [0.04 0.51 dx dx; ...
       0.51 0.51 dx dx; ...
       0.04 0.05 dx dx; ...
       0.51 0.05 dx dx];

figure(1); clf;
%set(gcf,'Position',[1203 388 905 954]);
set(gcf,'Position',[9 5 707 797]);


TM = PRTG.TM;
YRPLT = datevec(TM);
DV   = datevec(TM);
nrc  = length(TM);
TMplt = TM(1:dtt:end);  % dates to plot


nplt = size(YRPLT,1);
nxplt = length(TMplt);
XPLT = YRPLT(1:dtt:end,:);
iS = 1;
iE = nxplt;

if iStrt>0
  iS = min(find(XPLT(:,1)>=iStrt));
end
if iEnd>0
  iE = max(find(XPLT(:,1)==iEnd));
end

if s_fig>1 
  idnmb = datenum(s_fig,1,1);
  iT = min(find(TMplt>=idnmb));
  cc=iT;
  yr0=s_fig;

  while yr0==s_fig
%    fgnm = sprintf('%sarc08110_GrSh_prtlcs_%4.4i.png',pthfig,cc);
    fgnm = sprintf('%sarc08110_LbrGrSh_prtlcs_%4.4i.png',pthfig,cc);
    if ~exist(fgnm,'file'); break; end;
    cc=cc+1;
    dv0=datevec(TMplt(cc));
    yr0=dv0(1);
  end
  if yr0>s_fig
    fprintf('All frames exist for year %i\n',s_fig);
    return;
  end
  
  cc=cc-1;
  iS=cc+1;
  fprintf('Requested year: %i, Satart from %s, frame cc=%i\n\n',s_fig,datestr(TMplt(iS)),iS);
  
  
%keyboard
end

fprintf('Frames for %s - %s\n',datestr(TMplt(iS)),datestr(TMplt(iE)));

for iyr=iS:iE
  dnmb=TMplt(iyr);
  dv1=datevec(TMplt(iyr));
  yr0=dv1(1);

  itt=find(TM==dnmb);
  tic;	
  fprintf('Plotting %i/%i/%i\n',dv1(1:3));
  
%  id1=dnmb-dtpth;
%  it1=min(find(TM>=id1));
%  it2=itt;
    
  figure(1); clf;
  for ilr=1:4
    pos=POS(ilr,:);
    axes('Position',pos);
				pcolor(Hmsk); shading flat;
				colormap([0 0 0; 1 1 1]);
				freezeColors;
				caxis([0 1]);

				hold on;
				
				contour(HH,[-8000:1000:-100],'Color',[0.8 0.8 0.8]);
				contour(HH,[-800 -800],'Linewidth',1.6,'Color',[0.5 0.5 0.5]);
				
%		keyboard
				lr=SLR(ilr);
				IZ=find(PRTG.ZL==lr);
		% Plot particles
    xpz = PRTG.XP(IZ,itt);
    ypz = PRTG.YP(IZ,itt);
    II = inpolygon(xpz,ypz,IGR(:,1),IGR(:,2));
    nin = length(find(II==1));
				plot(xpz,ypz,'.',...
							'Markersize',13,...
							'Color',[0 0.4 0.8]);

%
% WLabr Shelf
    IZ=find(PRTL.ZL==lr);
  % Plot particles
    xpz = PRTL.XP(IZ,itt);
    ypz = PRTL.YP(IZ,itt);
    II = inpolygon(xpz,ypz,IGR(:,1),IGR(:,2));
    nin = length(find(II==1));
    plot(xpz,ypz,'.',...
       'Markersize',13,...
       'Color',[0.8 0.4 0]);

				plot(IGR(:,1),IGR(:,2),'-','Linewidth',1.5,...
									'Color',[0.6 0. 0.]);
				
				axis('equal');
				set(gca,'xlim',[xlim1 xlim2],...
					'ylim',[ylim1 ylim2]);
				set(gca,'xtick',[],'ytick',[]);
				
    clear stt		
    stt{1} = sprintf('%2.2i/%2.2i/%i',dv1(3),dv1(2),dv1(1));
    stt{2} = sprintf('Lr=%i, Z=%i m',SLR(ilr),ZLR(ilr));
    stt{3} = sprintf('In: %i prt',nin);
    text(660,929,stt,'Fontsize',10,'Color',[1 0.9 0]);
    
    
  end;

  bottom_text(txtb,'pwd',1,'Fontsize',7);
%    drawnow

  cc = cc+1;
    
  if s_fig>0
    fgnm = sprintf('%sarc08110_LbrGrSh_prtlcs_%4.4i',pthfig,cc);
    fprintf('Saving figure %s\n',fgnm);
    print('-dpng','-r150',fgnm);
  end
  %fprintf('1 record: %8.5f min\n\n',toc/60);
		     %    keyboard

  fprintf('1 record processed: %5.2f min\n',toc/60);

end







