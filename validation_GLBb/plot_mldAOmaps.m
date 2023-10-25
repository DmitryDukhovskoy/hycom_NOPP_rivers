% Plot MLD or ILD in the Arctic Ocean/N. Atl.
% 2D fields - spatial maps
% calculated based on Kara et al. (2000)
% or other approaches 
% MLD calculation is done in mld_ArcticOcean.m
%
% Use monthly mean fields
% Note that original global HYCOM GLBb0.08 outputs 
% have been "remapped" onto ARCc0.08 grid
% see ../REMAP_ARCc/remap_gridGLB2ARC_archv
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

s_fig  = 0; % =0 do not save figure; =1: save figure
%year = 1993;
%im=1;       % month
f_ild = 0;  % plot ILD maps

rg = 9806;  % convert pressure to depth, m
%dT = 0.2; % T change to calculate d(rho) for MLD - Kara et al., 2003 & 2000
dT = 0; % T change to calculate d(rho) for MLD - Kara et al., 2003 & 2000
          % dT=0 - uses my algorithm for MLD (see sub_mld_DD.m)

pthbin  = '/Net/mars/ddmitry/hycom/GLBb2ARCc0.08/mnth_mean/'; 
pthmat  = '/Net/mars/ddmitry/hycom/GLBb0.08/data_mld/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/'; 
pthfig  = '/Net/mars/ddmitry/hycom/GLBb0.08/fig_MLD/';

TV = '07';
ftopo = sprintf('%sdepth_ARCc0.08_%s.nc',pthtopo,TV); % 
HH  = nc_varget(ftopo,'Bathymetry'); 
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(HH); 


for yr=1993:1993
  for im=1:1
    year = yr;
    fprintf('Year %i, Month %i\n',yr,im);
%    fmat = sprintf('%sHYCOManls_monthlyMLD_dT%i%i',pthmat,dT*10,year);
    fmat = sprintf('%sHYCOManls_monthlyMLD_dT%2.2i_%i',pthmat,dT*10,year);
    if dT==0
      fmat = sprintf('%sHYCOManls_monthlyMLD_DD_%i',pthmat,year);
    end      
    fprintf('Loading %s\n',fmat);
    load(fmat);

    c1 = -1200;
    c2 = 0;
    
    if im>5  & im<10
      c1=-50;
      c2=0;
    end
    
    nint = 150;
    CMP = colormap_mld1(nint,c1,c2);
    cmp = CMP.colormap;
    cnt = CMP.intervals;

    mld=MILD(im).MLD;
    

    if f_ild>0 % plot Isothermal Layer Depth
      ild=MILD(im).ILD;
      figure(1); clf;
      pcolor(ild); shading flat;
      caxis([-1000 0]);
      hold on;
      contour(HH,[0 0],'k');
      axis('equal');
      set(gca,'xlim',[200 1400],...
	      'ylim',[10 1850]);
      clb=colorbar;
      %set(clb,'TickLength',0.047,...
      %	'Position',[0.82, 0.2, 0.025, 0.6],...
      %	'Fontsize',11);
      set(clb,'TickLength',0.045,...
	      'Position',[0.81, 0.15, 0.022, 0.7],...
	      'Fontsize',11);

      stt = sprintf('GLBb0.08, ILD, %i/%i',year,im);
      title(stt);
      btx='mld_ArcticOcean.m';
      bottom_text(btx,'pwd',1);
    end


    figure(2); clf;
    pcolor(mld); shading flat;
    colormap(cmp);
    caxis([c1 c2]);
    hold on;
    contour(mld,[-50 -50],'b');
    contour(mld,[-30 -30],'c');
    contour(mld,[-20 -20],'g');
    contour(HH,[0 0],'Color',[0.95 0.95 0.95]);
    axis('equal');
    set(gca,'xlim',[200 1400],...
	    'ylim',[10 1850],...
	    'xtick',[],...
	    'ytick',[],...
	    'Color',[0 0 0]);
    clr=[0.9 0.9 0.9];
    plot_gridlines(90,10,1,clr,LON,LAT);

    clb=colorbar;
    set(clb,'TickLength',0.04,...
	    'Position',[0.81, 0.15, 0.022, 0.7],...
	    'Fontsize',11);

    stt = sprintf('GLBb0.08, MLD, %i/%i',year,im);
    title(stt);
    btx='plot_mldAOmaps.m';
    bottom_text(btx,'pwd',1);

    if s_fig>0
      fgnm = sprintf('%sHYCOM_GLBb008_MLD_dT%2.2i_%i%2.2i',...
		     pthfig,round(dT*10),year,im);
      if dT == 0
        fgnm = sprintf('%sHYCOM_GLBb008_MLD_DD_%i%2.2i',...
		     pthfig,year,im);
      end	
      fprintf('Saving %s\n',fgnm);
      print('-dpng','-r200',fgnm);
    end
  end % month
end   % year






