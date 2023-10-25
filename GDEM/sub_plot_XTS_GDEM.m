function sub_plot_XTS_GDEM(ToPlot,Hgdem,zz,alat,elonM,pthfig,s_fig,Llon,Rlon,Fld,units,c1,c2);
 
% plots a cross section of ToPlot,at depths zz 
% plot is labeled as 'var' (variable name), with units 'units'
% contour levels c1, c2
%
% along left branch Llon, right branch Rlon, crossing the equator

fprintf('Plotting cross-section for %s ...\n',Fld);

lat0=65

% for the left branch, latitudes increase; 
L_lat=alat(alat>=lat0);
% for the right branch, latitudes decrease; 
R_lat=flipud(alat(alat>=lat0));
% string left and righr branches together
LR_lat=[L_lat;R_lat(2:end)];
% NB: Check if the pole is repeated, or we are on both sides of it?


% left and right branches for the variable being plotted
 L_ToPlot=squeeze(ToPlot(:,min(find(alat>=lat0)):end,find(elonM==Llon))); 
 R_ToPlot=fliplr(squeeze(ToPlot(:,min(find(alat>=lat0)):end,find(elonM==Rlon))));
% left and right branches for topography
 L_Hgdem=squeeze(Hgdem(min(find(alat>=lat0)):end,find(elonM==Llon)));
 R_Hgdem=flipud(squeeze(Hgdem(min(find(alat>=lat0)):end,find(elonM==Rlon))));
% string together left and right branches
 LR_ToPlot=[L_ToPlot,R_ToPlot(:,2:end)];
 LR_Hgdem=[L_Hgdem;R_Hgdem(2:end)];
% 


 for k=1:78
        LR_ToPlot(k,LR_Hgdem>=zz(k))=nan;
 end


nint= 200;
CMP = create_colormap2_1(nint,c1,c2);
cmp0= CMP.colormap;
nav = 10;
cmp = smooth_colormap(cmp0,nav);
cnt = CMP.intervals;

xx=[1:size(LR_lat)];
Depth=zz;


  nf=1;
  figure(nf); clf;
  

nint=7;
cnt=(c1:(c2-c1)/nint:c2);
colormap(jet(length(cnt)));

cmp=[];
nsint2=4;
jjs=clrmp_zebra(cnt,cmp,nsint2);


pcolor(xx,Depth,LR_ToPlot);
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

  %HT=title(['GLBb0.08, ',Fld,', ',int2str(yr1(decade)),'-',int2str(yr2(decade))]);
  %set(HT,'FontSize',12);

hold on;

  set(gca,'xtick',[100:100:800],'tickdir','out');
  set(gca,'ytick',[-4500:500:0],'tickdir','out');

  contour(xx,Depth,LR_ToPlot,[-0.75 -0.75],'k','LineWidth',1);
  contour(xx,Depth,LR_ToPlot,[0 0],'k','LineWidth',1);
  contour(xx,Depth,LR_ToPlot,[-1. -1.],'k','LineWidth',1);


  ctl=sprintf('GDEM4, %s(%s), cross-section at [%iE] and [%iE]',Fld,units,Llon,Rlon)
  title(ctl,'Fontsize',13);
  
  if s_fig>0
    fgnm=sprintf('%sGDEM4_xsect_%s_%i_%i',pthfig,Fld,Llon,Rlon);
    fprintf('Saving %s\n',fgnm);
    print('-dpng','-r200',fgnm);
  end

  
  
  



return
