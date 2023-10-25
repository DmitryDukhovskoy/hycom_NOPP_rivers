function sub_plot_bath3(HH,LON,LAT,fn,xyl);

xl1 = xyl(1,1);
yl1 = xyl(1,2);
xl2 = xyl(2,1);
yl2 = xyl(2,2);

figure(fn); clf;
set(gcf,'Position',[1400 345 951 986]);

%LMSK = HH*0;
%LMSK(HH<0)=1;
%lcmp = [0. 0. 0.; 1 1 1];
%pcolor(LMSK); shading flat;
%hold on;
%colormap(lcmp);
%freezeColors;

c1 = -5000;
c2 = 0;
CMP = create_colormap_bath(400,c1,c2);
cmp = flipud(CMP.colormap);

HH(HH>0)=nan;
pcolor(HH); shading flat;
hold on;
colormap(cmp);
caxis([c1 c2]);

%contour(HH,[0 0],'Color',[0.6 0.6 0.6]);
%contour(HH,[-500 -500],'Color',[0.9 0.7 0]);
%contour(HH,[-500 -500],'Color',[0.5 0.5 0.5]);
%contour(HH,[-8000:500:-100],'Color',[0.8 0.8 0.8]);
axis('equal');
set(gca,'xlim',[xl1 xl2],'ylim',[yl1 yl2]);
set(gca,'Color',[0 0 0]);

%set(gcf,'Position',[666 98 1507 1241]);

clb = colorbar;
set(clb,'Position',[0.85 0.1 0.02 0.82],...
        'Fontsize',14);

% Plot Lon/Lat;
dmm=LON;
dmm(dmm<-170)=nan;
dmm(dmm>170)=nan;
dmm(LAT>88) = nan;
%contour(dmm,[-180:30:180],'Color',[0.6 0.6 0.6]);
dmm=LON;
I=find(dmm<0);
dmm(I)=dmm(I)+360;
dmm(dmm>200)=nan;
dmm(dmm<160)=nan;
dmm(LAT>88)=nan;
%contour(dmm,[150:30:190],'Color',[0.6 0.6 0.6]);
%contour(LAT,[40:10:88],'Color',[0.6 0.6 0.6]);

return
