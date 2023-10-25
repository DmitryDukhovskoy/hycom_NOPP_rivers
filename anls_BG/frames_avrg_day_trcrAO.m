% Plot monthly-mean, depth-integrated 
% tracer conc fraction
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
pthfig  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/frames_day_trcrAO/',...
		  expt);
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
c2 = 3;  %4
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

POS(2,:) = [0.08 0.52 0.46 0.46];
POS(3,:) = [0.48 0.52 0.46 0.46];
POS(4,:) = [0.08 0.03 0.46 0.46];
POS(5,:) = [0.48 0.03 0.46 0.46];
    
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
    cc = cc+1;

    if cc<=fr_last,
      fprintf('Frame %4.4i need %4.4i, skipping ...\n',cc,fr_last);
      continue;
    end
    
    close all
    ff=figure('visible','off');
    
    for nTr = 2:5
      for ilv=1:1
        fprintf('Plotting: %2.2i/%2.2i/%2.2i, Tracer %i, Level=%i\n',...
	      dv(1:3),nTr,ilv);
        Tr = double(squeeze(TRCR(nTr).TR_lr1(idd,:,:)));
	Tr(Tr<=0)=nan;
	lTr=log(Tr);
	lTr(isnan(lTr))=-1000;
	lTr(HH>=0)=nan;

	if ilv==1
	  lvl='0-50m';
	else
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
	
        stl = sprintf('ARCc0.08-%3.3i, %2.2i/%2.2i/%2.2i, Trcr# %i, %s',...
		    expt,dv(1:3),nTr,lvl);
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
	
      end % ilv
    end  % nTr
% keyboard   
    if s_fig>0
      txtb='frames_avrg_day_trcrAO.m';
      bottom_text(txtb,'pwd',1,'Fontsize',8);
      fgnm = sprintf('%smnth_trcrLev%2.2i_%5.5i',pthfig,ilv,cc);
      fprintf('Saving figure %s\n',fgnm);
%      print('-djpeg','-r300',fgnm);
      print('-dpng','-r300',fgnm);
    end
    fprintf('1 record: %8.5f min\n\n',toc/60);

  end; % days
end   % year loop





