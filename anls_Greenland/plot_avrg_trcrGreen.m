% Plot monthly-mean, depth-integrated 
% tracer concentrations - kg/m3
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

s_fig = 1;
%dref = datenum(2006,1,15); % reference year/month, = 0 - no reference 
%ilev = 1; %=1 - mixed layer, 0-50; =2 - 50-150m

regn = 'ARCc0.08';
expt = 110;  
%pthfig  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/fig_trac/',expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);


ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

xlim1 = 20;
xlim2 = nn;
ylim1 = 5;
ylim2 = 2000;

hmsk=HH;
hmsk(HH<0)=nan;

nint = 360;
c1 = -1;
c2 = 1;
CMP = create_colormap2_3(nint,c1,c2);
cmp = CMP.colormap;
cmp = smooth_colormap(cmp,18);
cmp(1,:) = [1 1 1];
cnt = CMP.intervals;

for iyr=1993:2016
  pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_mnth%i/',...
		   regn,expt,iyr);
  if ~exist(pthfig,'dir')
    scm = sprintf('mkdir -pv %s',pthfig);
    system(scm);
  end
  

  for imo=7:7
    fmat = sprintf('%strcr_dpthav_%i%2.2i.mat',pthmat,iyr,imo);
    fprintf('Loading %s\n',fmat);
    load(fmat);
    
    for nTr = 1:5
      figure(1); clf;
%      ff=figure('visible','off');
      
      if nTr == 5 % Bering Str. - higher FW flux
	c1=-2;
	c2=3;
      else
	c1=-3;
	c2=2;
      end
      
      
      for ilv=1:2
	fprintf('Plotting: %i/%2.2i, Tracer %i, Lev %i\n',...
		iyr,imo,nTr,ilv);
        Tr = squeeze(TRCR(ilv).TR(nTr,:,:));
	Tr(Tr<=0)=nan;
	lTr=log(Tr);

	if ilv==1
	  pst = [0.03 0.08 0.4 0.85];
	  lvl='0-50m';
	else
	  pst = [0.47 0.08 0.4 0.85];
	  lvl='50-150m';
	end
	
	axes('Position',pst);
	pcolor(lTr); shading flat;
        hold on;
        contour(HH,[0 0],'k','Linewidth',1);
        contour(HH,[-5000:1000:-100],'Color',[0.7 0.7 0.7],'Linewidth',0.8);
        caxis([c1 c2]);
        colormap(cmp);
	axis('equal');
	set(gca,'xlim',[xlim1 xlim2],...
		'ylim',[ylim1 ylim2],...
		'Color',[0. 0. 0.]);
	set(gca,'xtick',[],'ytick',[]);
	clr=[0.9 0.9 0.9];
	plot_gridlines(45,10,0.5,clr,LON,LAT);
	stl = sprintf('%i/%2.2i, MAv, Trcr# %i, log(Conc), %s',...
		      iyr,imo,nTr,lvl);
	title(stl,'Fontsize',8);
	
	if ilv == 2
          clb=colorbar;
          set(clb,'TickLength',0.044,...
		'Position',[0.9, 0.2, 0.02, 0.6],...
		'Fontsize',10);
%          set(clb,'Ticks',[-1:0.2:1],...
%		'TickLength',0.045,...
%		'Position',[0.91, 0.2, 0.02, 0.6],...
%		'Fontsize',10);
        end	  
	
      end
      
      if s_fig>0
	txtb='plot_avrg_trcr.m';
	bottom_text(txtb,'pwd',1);
	fgnm = sprintf('%smnth_trcr%i_%i%2.2i',pthfig,nTr,iyr,imo);
	fprintf('Saving figure %s\n',fgnm);
	print('-dpng','-r250',fgnm);
      end
      
    end
    
  end % month
end   % year




