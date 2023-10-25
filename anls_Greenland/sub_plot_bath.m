function sub_plot_bath(HH,LON,LAT,fn,domname);
% Plot bath map for sections
switch(domname),
 case('NorthAtl')
  xl1=1;
  xl2=1100;
  yl1=10;
  yl2=1050;
 otherwise
  xl1=300;
  xl2=1300;
  yl1=50;
  yl2=1100;
end

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
contour(HH,[-500 -500],'Color',[0.5 0.5 0.5]);
contour(HH,[-8000:1000:-100],'Color',[0.9 0.9 0.9]);
axis('equal');
set(gca,'xlim',[xl1 xl2],'ylim',[yl1 yl2]);
%set(gcf,'Position',[666 98 1507 1241]);

% Plot Lon/Lat;
dmm=LON;
dmm(dmm<-170)=nan;
dmm(dmm>170)=nan;
contour(dmm,[-180:10:180],'Color',[0.6 0.6 0.6]);
dmm=LON;
I=find(dmm<0);
dmm(I)=dmm(I)+360;
dmm(dmm>190)=nan;
dmm(dmm<170)=nan;
contour(dmm,[170:10:190],'Color',[0.6 0.6 0.6]);
contour(LAT,[40:10:88],'Color',[0.6 0.6 0.6]);

return