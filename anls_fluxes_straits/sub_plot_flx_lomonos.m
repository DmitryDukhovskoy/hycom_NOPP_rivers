function sub_plot_flx_lomonos(HH,SCT,iFmn,iFmax,iFmin,fn);

hmsk=HH;
hmsk(HH<0)=nan;
cll = colormap([0.65 0.65 0.65]);

xl1 = 600;
xl2 = 1200;
yl1 = 1000;
yl2 = 1850;

% Plot Mean FWflux , upper 50m
figure(fn); clf;
pcolor(hmsk); shading flat;
colormap(cll);
freezeColors;
%contour(HH,[0 0],'k');
hold on;

contour(HH,[-5000:500:-10],'Color',[0.75 0.75 0.75]);
Ept=SCT.End_Pnts;
for ig=1:2
  plot([Ept(ig,1) Ept(ig+1,1)],[Ept(ig,2) Ept(ig+1,2)],...
     '-','Color',[0.2 0.2 0.2],...
     'Linewidth',1.6);
end

plot(iFmax(:,1),iFmax(:,2),'-','Linewidth',1.8,'Color',[0 0.7 0.8]);
plot(iFmin(:,1),iFmin(:,2),'-','Linewidth',1.8,'Color',[0 0.7 0.8]);
plot(iFmn(:,1),iFmn(:,2),'r-','Linewidth',2);
  
axis('equal');
set(gca,'xlim',[xl1 xl2],...
	'ylim',[yl1 yl2],...
	'xtick',[],...
	'ytick',[]);

return