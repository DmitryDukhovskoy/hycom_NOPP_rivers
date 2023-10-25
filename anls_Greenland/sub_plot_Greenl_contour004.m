function sub_plot_Greenl_contour004(HH,LON,LAT,fn,GC);
% Plot contour for flux analysis
Hg = HH;

Hg(1:760,:)=nan;
Hg(2200:end,:)=nan;
Hg(:,2100:end)=nan;
Hg(:,1:900)=nan;

fprintf('Plotting Greenland section/contour ...\n');
fn=2;
domname = '0';
sub_plot_bath(Hg,LON,LAT,fn,domname);
contour(Hg,[-100:10:-5],'Color',[0.85 0.85 0.85]);
contour(Hg,[-3500:100:-10],'Color',[0.6 0.6 0.6]);
contour(Hg,[-4000:500:-50],'Color',[0.25 0.25 0.25]);
IIs = GC.cntr_Iindx;
JJs = GC.cntr_Jindx;
x   = GC.Distance_m*1e-3; % m->km

plot(IIs,JJs,'b-','Linewidth',2);
for km=0:500:max(x)
  d=abs(x-km);
  i0=find(d==min(d));
  if km==0
     plot(IIs(i0),JJs(i0),'r.','Markersize',14);
     plot(IIs(i0),JJs(i0),'rd','Markersize',6);
  else
    plot(IIs(i0),JJs(i0),'r.','Markersize',11);
  end
  text(IIs(i0),JJs(i0),sprintf('%6.1f km',km),'Fontsize',14);
end

set(gca,'xlim',[900 2100],...
	'ylim',[760 2200],...
	'xtick',[],...
	'ytick',[]);

title('Greenl. Contour (~800m) for shelf-basin flux analysis');


return
