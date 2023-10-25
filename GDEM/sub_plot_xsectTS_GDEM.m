function sub_plot_xsectTS_GDEM(ToPlot,Depth,Fld,IJ,decade,yr1,yr2,c1,c2,pthfig,id)
%
% ToPlot: variable whose cross-section this is
% Depth: depths
% Field name of variable with units
% yr1 and yr2: first and last year
% IJ number of points from left to right
% c1, c2: contour levels

  nf=1;
  figure(nf); clf;

nint=7;
cnt=(c1:(c2-c1)/nint:c2);
colormap(jet(length(cnt)));

cmp=[];
nsint2=4;
jjs=clrmp_zebra(cnt,cmp,nsint2);


xx=(1:length(IJ));
pcolor(xx,Depth,ToPlot);
shading interp;
caxis([c1 c2]);

colormap(jjs);

set(gca,'ylim',[-4700 0]);
set(gca,'Color',[0 0 0],'tickdir','out');
set(gcf,'InvertHardcopy','off')

sfig=1

    hght=[];
    lngth=[];
    mint=1;
    mbx=1;
    fsz=10;
    bxc='w';
    posc=[0.23 0.04 0.6 0.05];
    aend=1;
[az,axc]  = pcolorbar_zebra_horiz(jjs,cnt,hght,lngth,mint,fsz,bxc,posc,nsint2);
%  [az,axc]  = colorbar_horiz (jjs,cnt,hght,lngth,mint,fsz,bxc,posc,mbx,aend);

  HT=title(['GLBb0.08, ',Fld,', ',int2str(yr1(decade)),'-',int2str(yr2(decade))]);
  set(HT,'FontSize',12);

hold on;

set(gca,'xtick',[100:100:800],'tickdir','out');
set(gca,'ytick',[-4500:500:0],'tickdir','out');

contour(xx,Depth,ToPlot,[-0.75 -0.75],'k','LineWidth',1);
contour(xx,Depth,ToPlot,[0 0],'k','LineWidth',1);
contour(xx,Depth,ToPlot,[-1. -1.],'k','LineWidth',1);

  if sfig==1
    disp('Saving figure ...');
    fgnm=sprintf('%sGLBb008_xsect_%s_%s_%i-%i.png',pthfig,id,Fld,yr1(decade),yr2(decade));
    print('-dpng','-r200',fgnm)
  end;
