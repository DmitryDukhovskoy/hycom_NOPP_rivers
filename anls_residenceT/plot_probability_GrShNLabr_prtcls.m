% Plot trajectories of all partilces
% released on the W Gr Shelf and NW Labr Shelf
% by layers
% all output have to be preprocessed in
% plot_WLabrShprt_paths.m for all years
% plot_GrShprt_paths.m
% 
% The code is not finished for combining GrSh and Labr Sh particles
% Use it to plot Labr Sh probabilities of partilces
% for Gr Shelf can use this one or plot_probability_prtcls.m
% should be same codes
%
addpath /Net/Movies0/ddmitry/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /home/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers;
%addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
%addpath /Net/Movies0/ddmitry/MyMatlab/seawater
startup

close all
clear

f_plt=3; % =2 - plot Lab Shelf particles only
         % =3 - plot combined probability for Gr Sh and Labr Sh.

SLR=[10,15,23,31]; % depth levels of the particels
ZLR=[50,90,150,450];  % nominal depths of particles in deep ocean
nlr=length(SLR);

s_fig = 0;

% blkdat.input: target densities for the layers:
% upper 14 layers - fixed z-layers over the deep ocean
% 26.00   'sigma ' = layer 10 isopycnal target density (sigma units) - lr 5
% 30.65   'sigma ' = layer  A isopycnal target density (sigma units) - lr 15
% 35.20   'sigma ' = layer  I isopycnal target density (sigma units) - lr 23
% 36.70   'sigma ' = layer 22 isopycnal target density (sigma units) - lr 31

regn = 'ARCc0.08';
expt = 110;  
%pthfig  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/fig_trac/',expt);
pthfig  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/%s/%3.3i/fig_GrGprtcl/',regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/110/data_GG_prt/';
pthmat2 = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/110/data_WLabr_prt/';

txtb = 'plot_probability_GrShNLabr_prtcls.m';

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



fnmout=sprintf('%sGrShprt_comb.mat',pthmat);
fnmout2=sprintf('%sWLabrShprt_comb.mat',pthmat2);

% Combine all experiments with Lagrangian particles together 
% GrShelf and Labr shelf separately
% for missing dates - use the previous number and locations of particles
PPG = sub_GrSh(fnmout);  % Gr shelf
PPL = sub_GrSh(fnmout2); % Labr shelf

% Combine Gr Shelf and Labr Shelf:
% not finished combining all trajectories and time TM
if f_plt==3
		PPGL=[];
		xp1=PPG(1).XPcmb;
		yp1=PPG(1).YPcmb;
		zl1=PPG(1).ZLcmb;
		xp2=PPL(1).XPcmb;
		yp2=PPL(1).YPcmb;
		zl2=PPL(1).ZLcmb;
  PPGL(1).XPcmb=[xp1;xp2];
  PPGL(1).YPcmb=[yp1;yp2];
  PPGL(1).ZLcmb=[zl1;zl2];

  PRGL = sub_prtcl_stat(PPGL,HH,nlr,SLR,xlim1,xlim2,ylim1,ylim2);  % Labr Sh prob
  PR = PRGL;
end

% Statistics separatle for Labr and Gr shelves
%PRG = sub_prtcl_stat(PPG,HH,nlr,SLR,xlim1,xlim2,ylim1,ylim2); % Gr Sh probab
if f_plt==2
  PRL = sub_prtcl_stat(PPL,HH,nlr,SLR,xlim1,xlim2,ylim1,ylim2);  % Labr Sh prob
  PR = PRL;
end
%

% Time size should match - # of columns
PPGL(1).

Hmsk=HH*0;
Hmsk(HH<0)=1;

% Region of interest - double check
% with dS_FWC_timeseries_SubpolarGyre.m
IGR = [  430         729
         461         671
         559         660
         774         676
         781         608
         839         533
         855         530
         930         501
        1031         445
        1099         410
        1065         350
        1087         288
        1106         204
        1118         148
        1080         145
         444         145
         416         246
         367         492
         374         591
         379         729
        430       729];


  
CMP = create_colormap_WBYR(400,0,1);
cmp = CMP.colormap;

%
for ilr=1:nlr
  figure(ilr); clf;
  pcolor(Hmsk); shading flat;
  colormap([0 0 0; 1 1 1]);
  hold on;
  
  caxis([0 1]);
  freezeColors;

  dmm = PR(ilr).SMM;
  dmm(HH>=0)=nan;
  I0=find(dmm==0);
  dmm(I0)=nan;
  dmm=log(dmm);
  dmm(I0)=-80;
  
  pcolor(dmm); shading flat;
		colormap(cmp);
  caxis([-18 -11]);

  plot(IGR(:,1),IGR(:,2),'-','Linewidth',1.5,...
       'Color',[0.6 0. 0.]);
  contour(HH,[-1000 -1000],'Color',[0.6 0.6 0.6]);
  contour(HH,[-500 -500],'Linewidth',1.6,'Color',[0.4 0.4 0.4]);
  
  axis('equal');
  set(gca,'xlim',[xlim1 xlim2],...
	  'ylim',[ylim1 ylim2]);
  set(gca,'xtick',[],'ytick',[]);
    
  clb = colorbar;
  set(clb,'Position',[0.91 0.11 0.017 0.81],...
          'TickLength',0.025,...
          'Fontsize',14);

  if f_plt==2
    stt = sprintf('Prob Lagr prtcl in grid cell, log2, %i m',ZLR(ilr));
  elseif f_plt==3
    stt = sprintf('Prob GrSh & NLabrSh prtcl in grid cell, log2, %i m',ZLR(ilr));
  end
  title(stt);
    

  bottom_text(txtb,'pwd',1);
    
%  if s_fig>0
%    fgnm = sprintf('%s008-110_NLabrSh_prtcl_lr%i',pthfig,ilr);
%    fprintf('Saving figure %s\n',fgnm);
%    print('-dpng','-r200',fgnm);
%  end
  %fprintf('1 record: %8.5f min\n\n',toc/60);
		     %    keyboard

end





