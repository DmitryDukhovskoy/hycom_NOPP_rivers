% Plot monthly-mean, depth-integrated 
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

s_fig = 1;
cc=0; % last SAVED frame - if need to restart from frame # 
      % the next frame will be cc+1
      % = 0 - start from beginning

%dref = datenum(2006,1,15); % reference year/month, = 0 - no reference 
%ilev = 1; %=1 - mixed layer, 0-50; =2 - 50-150m

regn = 'ARCc0.08';
expt = 110;  
pthfig  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/frames_trcr/',expt);
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

nint = 320;
c1 = -3; %-4
c2 = 2;  %4
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
for iyr = 1993:2015
  for imo=1:12
    ck=ck+1;
    YRPLT(ck,1)=iyr;
    YRPLT(ck,2)=imo;
  end
end
nrc=ck;

%cc=0;
for ic=cc+1:nrc
  iyr = YRPLT(ic,1);
  imo = YRPLT(ic,2);
  
%  pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_mnth%i/',...
%		   regn,expt,iyr);
%  if ~exist(pthfig,'dir')
%    scm = sprintf('mkdir -pv %s',pthfig);
%    system(scm);
%  end
    cc=cc+1;
    fmat = sprintf('%strcr_dpthav_%i%2.2i.mat',pthmat,iyr,imo);
    if ~exist(fmat,'file')
      cotinue
    end
    
%    if iyr~=2008 | imo~=7, continue; end;
tic;    
    fprintf('cc=%i, Loading %s\n',cc,fmat);
    load(fmat);
    
    close all
    ff=figure('visible','off');
%    figure(1); clf;
    POS(2,:) = [0.08 0.52 0.45 0.45];
    POS(3,:) = [0.48 0.52 0.45 0.45];
    POS(4,:) = [0.08 0.03 0.45 0.45];
    POS(5,:) = [0.48 0.03 0.45 0.45];
    
    for nTr = 2:5
%      figure(1); clf;
      
%      if nTr == 5 % Bering Str. - higher FW flux
%	c1=-4;
%	c2=4;
%      else
%	c1=-4;
%	c2=4;
%      end
      for ilv=1:1
	fprintf('Plotting: %i/%2.2i, Tracer %i, Lev %i\n',...
		iyr,imo,nTr,ilv);
        Tr = squeeze(TRCR(ilv).TR(nTr,:,:));
	Tr(Tr<=0)=nan;
	lTr=log(Tr);
	lTr(isnan(lTr))=-1000;
	lTr(HH>=0)=nan;

	if ilv==1
%	  pst = [0.03 0.08 0.4 0.85];
	  lvl='0-50m';
	else
%	  pst = [0.47 0.08 0.4 0.85];
	  lvl='50-150m';
	end
	
        pst = POS(nTr,:);
	axes('Position',pst);
	pcolor(lTr); shading flat;
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
	plot_gridlines(45,20,0.5,clr,LON,LAT);
	stl = sprintf('%i/%2.2i, Trcr# %i, %s',...
		      iyr,imo,nTr,lvl);
	title(stl,'Fontsize',8);
	
	if nTr == 3 
          clb=colorbar;
          set(clb,'TickLength',0.02,...
		'Position',[0.91, 0.55, 0.02, 0.4],...
		'Fontsize',10);
        end	  
	if nTr == 5 
          clb=colorbar;
          set(clb,'TickLength',0.02,...
		'Position',[0.91, 0.08, 0.02, 0.4],...
		'Fontsize',10);
        end	  
      end
      
    end
% keyboard   
    if s_fig>0
      txtb='frames_avrg_trcr.m';
      bottom_text(txtb,'pwd',1);
      fgnm = sprintf('%smnth_trcrLev%2.2i_%5.5i',pthfig,ilv,cc);
      fprintf('Saving figure %s\n',fgnm);
%      print('-djpeg','-r300',fgnm);
      print('-dpng','-r300',fgnm);
    end
    fprintf('1 record: %8.5f min\n\n',toc/60);

    
end   % time loop




