function sub_plot_bath2(HH,LON,LAT,fn,xyl);

xl1 = xyl(1,1);
yl1 = xyl(1,2);
xl2 = xyl(2,1);
yl2 = xyl(2,2);

figure(fn); clf;
LMSK = HH*0;
LMSK(HH<0)=1;
lcmp = [0.4 0.4 0.4; 1 1 1];
pcolor(LMSK); shading flat;
colormap(lcmp);
freezeColors;

hold on;
%contour(HH,[0 0],'Color',[0.6 0.6 0.6]);
%contour(HH,[-500 -500],'Color',[0.9 0.7 0]);
%contour(HH,[-500 -500],'Color',[0.5 0.5 0.5]);
%contour(HH,[-8000:1000:-100],'Color',[0.9 0.9 0.9]);
axis('equal');
set(gca,'xlim',[xl1 xl2],'ylim',[yl1 yl2]);
%set(gcf,'Position',[666 98 1507 1241]);

% Plot Lon/Lat;
dmm=LON;
dmm(dmm<-170)=nan;
dmm(dmm>170)=nan;
dmm(LAT>88) = nan;
contour(dmm,[-180:20:180],'Color',[0.6 0.6 0.6]);
dmm=LON;
I=find(dmm<0);
dmm(I)=dmm(I)+360;
dmm(dmm>200)=nan;
dmm(dmm<160)=nan;
dmm(LAT>88)=nan;
contour(dmm,[170:20:190],'Color',[0.6 0.6 0.6]);
contour(LAT,[40:10:88],'Color',[0.6 0.6 0.6]);

return