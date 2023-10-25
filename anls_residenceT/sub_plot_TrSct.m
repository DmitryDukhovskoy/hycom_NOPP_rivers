% Plot vertical section
% Tracer and U contours
% on Greenland shelves
function sub_plot_TrSct(isc,DX,ZZ,lmVTr,xBtm,yBtm,xU,mUU,c1,c2,yl1,xl1,xl2,CTR);

% log scale, m3 of GFWA in 1 m3 of sea water, i.e. liters of GFWA
if isempty(c1)
  c1=-8;
  c2=-4.5;
end
nint=300;
CMP=create_colormap8(nint,c1,c2);
cmp=CMP.colormap;
cnt=CMP.intervals;


figure(isc); clf;
set(gcf,'Position',[1156 602 957 698]);
axes('Position',[0.09 0.32 0.86 0.6]);
hold on;
pcolor(DX,ZZ,lmVTr); shading interp

contour(xU,ZZ,mUU,[0 0],'k-','Linewidth',1.6);
if isempty(CTR)
  ctr1=[0.02:0.02:0.4];
  ctr2=[-0.4:0.02:-0.001];
else 
  ctr1=CTR.Upos;
  ctr2=CTR.Uneg;
end
contour(xU,ZZ,mUU,ctr1,'k-','Color',[0.2 0.2 0.2]);
contour(xU,ZZ,mUU,ctr2,'k--','Color',[0.2 0.2 0.2]);

colormap(cmp);
caxis([c1 c2]);

patch(xBtm,yBtm,[0 0 0]);

clb=colorbar('SouthOutside');
set(clb,'Position',[0.09 0.2 0.86 0.03],...
								'TickLength',0.03);

set(gca,'tickdir','out',...
								'xlim',[xl1 xl2],...
								'ylim',[yl1 0],...
								'Fontsize',14);



return
