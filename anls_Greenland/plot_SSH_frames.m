% Plot ssh near Denmark strait
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers
startup;

close all
clear

s_fig = 1; % =0 - do not save mat file


rg = 9806;
yr1=2005;
yr2=2005;

fprintf('Years to extract: %i-%i\n',yr1,yr2);

ilim1 = 600;
ilim2 = 910;
jlim1 = 330;
jlim2 = 660;

% Experiments:
% 102 - test simulation, HYCOM GLBb0.08 nest: 1 file kept constant (1993,1,1)

regn = 'ARCc0.08';
expt = 110;
pthfig  = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/fig_ssh_frames/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
fmat    = sprintf('%sLapl_SSH_BG.mat',pthmat);


ftopo = sprintf('%sdepth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

[DX,DY]=sub_dx_dy(LON,LAT);
[II,JJ]=meshgrid([1:nn],[1:mm]);
dmm = inpolygon(II,JJ,[ilim1, ilim2, ilim2, ilim1],[jlim1, jlim1, jlim2, jlim2]);
Irg = find(dmm==1 & HH<-200.);


hmsk=HH;
hmsk(HH<0)=nan;


cl1 = colormap_red(200);
cl2 = flipud(colormap_blue(200));
for i=1:5
  cl1(1,:)=[1 1 1];
  cl2(end-i+1,:)=[1 1 1];
end

clrmp = [cl2;cl1];
clrmp = smooth_colormap(clrmp,9);


cc=0;
dday=2;
figure(1);
set(gcf,'Position',[1675, 620, 854, 707]);
btx = 'plot_SSH_frames.m';

yr=yr1;
for iday=1:dday:365
	pthbin = sprintf('/nexsan/archive/%s_%3.3i/data/%i/',regn,expt,yr);

	dnmb=datenum(yr,1,1)+iday-1;
	DV=datevec(dnmb);
	imo=DV(2);
	fprintf('Processing %i/%2.2i/%2.2i\n',DV(1:3));

	fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
	finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);

  cc=cc+1;

	tic;
	[F,nn,mm,ll] = read_hycom(fina,finb,'srfhgt');
	F(F>1e6)=nan;
	E=squeeze(F/9.806);

  ssh = E-nanmean(E(Irg));
  c1=-0.25;
  c2=0.25;
 
  clf; 
  axes('Position',[0.1 0.1 0.8 0.8]);
  hold on;
  pcolor(ssh); shading('flat');
  colormap(clrmp);
  contour(HH,[-8000:1000:0],'color',[0.5, 0.5, 0.5]);
  contour(ssh,[0:0.05:0.8],'color',[0.9 0.9 0.9]);
  contour(ssh,[-0.8:0.05:-0.01],'--','color',[0.9 0.9 0.9]);

  caxis([c1 c2]);
  axis('equal');
  set(gca,'xlim',[ilim1 ilim2],...
          'ylim',[jlim1 jlim2],...
          'color',[0 0 0],...
          'xticklabel',[],...
          'yticklabel',[]);
  
  ctl = sprintf('0.08HYCOM-CICE-110, SSH, %2.2i/%2.2i/%4.4i',DV(3),DV(2),DV(1));
  title(ctl);

  hb = colorbar;
  set(hb,'TickLength',0.03,...
       'Position',[0.9 0.1 0.02 0.8],...
       'Fontsize',14);

  bottom_text(btx,'pwd',1);

  ffg = sprintf('%s008hycom_110_GrSSH_%3.3i.png',pthfig,cc);
  if s_fig>0
    fprintf('Saving %s\n\n',ffg);
    print('-dpng','-r200',ffg);
  end
  
  sprintf('Processing %6.2f min\n',toc/60);

end

