% Plot ratio of Eurasian tracers 
% monthly-mean, depth-integrated 
% tracer concentrations 
% The Data are extracted in extr_trcr_mnth.m
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

rg  = 9806;
hgg = 1e20; % 


ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

xlim1 = 660;
xlim2 = 1350;
ylim1 = 1380;
ylim2 = 1950;

hmsk=HH;
hmsk(HH<0)=nan;

nint = 320;
c1 = -8; %-4
c2 = 8;  %4
%CMP = create_colormap2_3(nint,c1,c2);
CMP = colormap_sclr4(nint,c1,c2);
cmp = CMP.colormap;
cnt = CMP.intervals;

ck = 0;
for iyr = 2013:2014
  for imo=9:9
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
%    ff=figure('visible','off');
    figure(1); clf;

  %
  % Read S
    expt='110';
    pthbin = sprintf('/nexsan/archive/ARCc0.08_%s/data/%4.4i/',expt,iyr);
    
    iday = datenum(iyr,imo,15)-datenum(iyr,1,1)+1;
    fina = sprintf('%s%s_archm.%4.4i_%3.3i_12.a',pthbin,expt,iyr,iday);
    finb = sprintf('%s%s_archm.%4.4i_%3.3i_12.b',pthbin,expt,iyr,iday);


    tic;
    [F,n,m,nlr] = read_hycom(fina,finb,'salin','r_layer',1);
    F(F>hgg)=nan;
    S=squeeze(F);


    ilv=1;
    ctr=0;
    for nTr = 3:4 % tr3 - Lena, Tr4 - Enisei
      fprintf('Plotting: %i/%2.2i, Tracer %i, Lev %i\n',...
	      iyr,imo,nTr,ilv);
      Tr = squeeze(TRCR(ilv).TR(nTr,:,:));
      ctr=ctr+1;
      if ctr==1
	Tr1 = Tr;
      else
	Tr2 = Tr;
      end
    end;
%    keyboard   
    Tr2(Tr2<1e-5)=nan;
    Tr1(Tr1<1e-5)=nan;
    rTr = log(Tr1./Tr2);
    
    if ilv==1
%	  pst = [0.03 0.08 0.4 0.85];
      lvl='0-50m';
    else
%	  pst = [0.47 0.08 0.4 0.85];
      lvl='50-150m';
    end

    pst = [0.08 0.08 0.82 0.85];
    axes('Position',pst);
    pcolor(rTr); shading flat;
    hold on;
%    contour(HH,[0 0],'k','Linewidth',1);
%    contour(HH,[-1000 -1000],'Color',[0.7 0.7 0.7],'Linewidth',0.8);
    contour(HH,[-5000:1000: -10],'Color',[0.7 0.7 0.7],'Linewidth',0.8);
    caxis([c1 c2]);
    colormap(cmp);
    
    contour(S,[30 30],'k','linewidth',1.6);
    contour(S,[10:1:34],'Color',[0.2 0.2 0.2]);
    axis('equal');
    set(gca,'xlim',[xlim1 xlim2],...
	    'ylim',[ylim1 ylim2],...
	    'Color',[0. 0. 0.]);
    set(gca,'xtick',[],'ytick',[]);
    clr=[0.9 0.9 0.9];
    plot_gridlines(30,10,0.5,clr,LON,LAT);
    stl = sprintf('%i/%2.2i, log(Tr3/Tr4), %s',...
		  iyr,imo,lvl);
    title(stl,'Fontsize',12);

    clb=colorbar;
    set(clb,'TickLength',0.02,...
	    'Position',[0.91, 0.2, 0.02, 0.7],...
	    'Fontsize',12);

      txtb='tracerEuras_ratio.m';
      bottom_text(txtb,'pwd',1);
% keyboard   
    if s_fig>0
%      fgnm = sprintf('%smnth_ratio_Tr3Tr4_%3.3i',pthfig,cc);
      fgnm = sprintf('%smnth_ratio_Tr3Tr4_%i%2.2i',pthfig,iyr,imo);
      fprintf('Saving figure %s\n',fgnm);
%      print('-djpeg','-r300',fgnm);
      print('-dpng','-r300',fgnm);
    end
    fprintf('1 record: %8.5f min\n\n',toc/60);

    
end   % time loop




