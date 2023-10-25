% Plot daily-mean, depth-averaged 
% tracer concentrations 
% The Data are extracted in anls_BG/extr_trcr_day.m
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers
startup;

close all
clear

ilv = 1; % plot level 1 (0-50m)
s_fig = 1;

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
%dref = datenum(2006,1,15); % reference year/month, = 0 - no reference 
%ilev = 1; %=1 - mixed layer, 0-50; =2 - 50-150m

regn = 'ARCc0.08';
expt = 110;  
%pthfig=sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/frames_day_trcrGreen/',...
%		  expt);
pthfig='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/110/figs_trcr1/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);

ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

% For Greenland:
%xlim1 = 50;
%xlim2 = nn;
%ylim1 = 5;
%ylim2 = 1500;

xlim1 = 50;
xlim2 = nn-200;
ylim1 = 5;
ylim2 = 1150;

hmsk=HH;
hmsk(HH<0)=nan;

nint = 320;
c1 = -3;
c2 = 2;
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

%cc=2;
ic1=1;
cc=0; % last saved frame - if need to restart from frame # 
        % = 0 - start from beginning
for ic=ic1:nrc
  iyr = YRPLT(ic);
  fmat = sprintf('%strcr_dpthav_daily_%4.4i.mat',pthmat,iyr);
  
  if ~exist(fmat,'file')
    fprintf('File is missing %s\n',fmat);
    continue
  end

%  cc=cc+1;
  fprintf('Loading %s\n',fmat);
  load(fmat);

% keyboard
  
tic;    
  POS(1,:) = [0.08 0.08 0.85 0.85];
  POS(2,:) = [0.06 0.52 0.45 0.45];
  POS(3,:) = [0.48 0.52 0.45 0.45];
  POS(4,:) = [0.06 0.03 0.45 0.45];
  POS(5,:) = [0.48 0.03 0.45 0.45];

  TM = TRCR(1).TM;
  ndy = length(TM);
%  for idd = 1:ndy
  for idd = 30:30
    dnmb = TM(idd);
    dv = datevec(dnmb);
    cc = cc+1;
    
    if cc<=fr_last,
      fprintf('Frame %4.4i need %4.4i, skipping ...\n',cc,fr_last);
      continue;
    end
    
    
    close all
%    ff=figure('visible','off');
    figure(1); clf;

    for nTr = 1:1 % nTr = 1 - Greenland

      if ilv~=1,
	fprintf('Change code for layer>1');
	error(' current code of ilr=1')
      end

      fprintf('Plotting: %2.2i/%2.2i/%2.2i, Tracer %i, Level=%i\n',...
	      dv(1:3),nTr,ilv);
      Tr = double(squeeze(TRCR(nTr).TR_lr1(idd,:,:)));
      Tr(Tr<=0)=nan;
      lTr=log(Tr);
      lTr(isnan(lTr))=-1000;
      lTr(HH>=0)=nan;

      if ilv==1
%	  pst = [0.03 0.08 0.4 0.85];
	lvl='0-50m';
      else
%	  pst = [0.47 0.08 0.4 0.85];
	lvl='50-150m'; %
      end
%keyboard
      pst = POS(nTr,:);
      axes('Position',pst);
      pcolor(lTr); shading flat;
      hold on;
      contour(HH,[0 0],'k','Linewidth',1);
      contour(HH,[-5000:1000:-1000],'Color',[0.3 0.3 0.3],'Linewidth',0.8);
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
      title(stl);

      clb=colorbar;
      set(clb,'TickLength',0.02,...
	    'Position',[0.885, 0.1, 0.02, 0.8],...
	    'Fontsize',13);
%keyboard	  
      if s_fig>0
	txtb='frames_avrg_day_trcrGreen.m';
	bottom_text(txtb,'pwd',1,'Fontsize',6);
	fgnm = sprintf('%sdayGrtrcrLev%2.2i_%5.5i',pthfig,ilv,cc);
	fprintf('Saving figure %s\n',fgnm);
	print('-dpng','-r250',fgnm);
    %      print('-djpeg','-r350',fgnm);
      end
      fprintf('1 record: %8.5f min\n\n',toc/60);

    end  % nTr

  end  % day
  
end   % time loop




