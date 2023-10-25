function sub_plot_ssh(E,HH,LN,LT,nf);
%
%
btx = 'sub_plot_ssh.m';

xl1 = 180;
xl2 = 1360;
yl1 = 400;
yl2 = 1800;

fvs = (nf<0);
nf = abs(nf);

if fvs
  figure(nf); clf;
  set(gcf,'Visible','off');
else
  figure(nf); clf;
end

Hmsk = HH*0;
Hmsk(HH>0)=1;
cmm = [0.7 0.7 0.7; 0.,0.,0.];

nint = 200;
c1 = -0.5;
c2 = 0.5;
cl2 = colormap_hot(100);
cl1 = colormap_cold(100);
cl1 = flipud(cl1);
cmp = [cl1;cl2];

axes('Position',[0.05 0.1 0.8 0.8]);
pcolor(Hmsk); shading flat;
colormap(cmm);
caxis([0 1]);
hold on;
freezeColors;

pcolor(E); shading flat;
hold on
colormap(cmp);
caxis([c1 c2]);

contour(E,[0:0.1:1],'Color',[0.3 0.3 0.3]);
contour(E,[-1:0.1:-0.01],'k--','Color',[0.3 0.3 0.3]);

dmm=LN;
LN(LN>140) = nan;
LN(LN<-140)= nan;
LN(LT>89)  = nan;
contour(LN,[-140:20:140],'color',[0.8 0.8 0.8],'linewidth',1);
LN = dmm;
LN(LN<0) = LN(LN<0)+360;
LN(LT>89)  = nan;
LN(LN>360) = nan;
LN(LN<140) = nan;
contour(LN,[160:20:200],'color',[0.8 0.8 0.8],'linewidth',1);
contour(LT,[40:10:89],'color',[0.8 0.8 0.8],'linewidth',1);
axis('equal');
set(gca,'xlim',[xl1 xl2],...
	'ylim',[yl1 yl2],...
	'xtick',[],...
	'ytick',[]);

ch = colorbar;
set(ch,'position',[0.82 0.1 0.018 0.8],...
       'TickLength',0.022);


return