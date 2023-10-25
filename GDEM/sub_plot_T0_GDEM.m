function sub_plot_TS_GDEM(ToPlot,pthfig,s_fig,HH,LON,LAT,season,when,c1,c2);
% YR - array of year to plot
%      NxM, if M=1 - no long-term averaging
%           if M=2, for each n=1:N, ADPTH is averaged
%             for years yr1- yr2
fprintf('Plotting T_SFC for %s,...',when(season,:));


target=0;
var=sprintf('T_SFC for %s',when(season,:));

nint= 200;
CMP = create_colormap2_1(nint,c1,c2);
cmp0= CMP.colormap;
nav = 10;
cmp = smooth_colormap(cmp0,nav);
cnt = CMP.intervals;



[mm,nn]=size(HH);
Lm=HH;
Lm(Lm<target)=nan;
Lm(Lm>=target)=1;
cbb=[0 0 0; 1 1 1];

  nf=1;
  figure(nf); clf;
  pcolor(Lm); shading flat;
  colormap(cbb);
  hold on;
  freezeColors;
  
  pcolor(ToPlot); shading flat;
  caxis([c1 c2]);
  colormap(cmp);
  
  contour(HH,[target target],'Linewidth',1,...
	  'Color',[0.9 0.9 0.9]);
  
  axis('equal');
  set(gca,'Color',[0.3 0.3 0.3]);
  set(gca,'xlim',[250 nn],...
	'ylim',[400 1950],...
	'xtick',[],...
	'ytick',[]);
  
  dlmb=20;
  dphi=10;
  clr=[0.9 0.9 0.9];
  plot_gridlines(dlmb,dphi, nf, clr, LON, LAT)
  
  hght=[];
  lngth=[];
  mint=20;
  mbx=mint;
  fsz=13;
  bxc='k';
  posc=[0.8 0.11 0.8 0.08];
  aend=0;
  [az,axc]  = colorbar_vert (cmp,cnt,hght,lngth,mint,fsz,bxc,posc,mbx,aend);

  ctl=sprintf('GDEM4, T SFC (C) for %s',when(season,:));
  title(ctl,'Fontsize',13);
  
  if s_fig>0
    fgnm=sprintf('%sGDEM4_T_SFC_%s',pthfig,when(season,:));
    fprintf('Saving %s\n',fgnm);
    print('-dpng','-r200',fgnm);
  end

  
  
  



return
