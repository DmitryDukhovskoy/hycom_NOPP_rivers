function sub_plot_regions(nf,FWB,HH,LON,LAT,nrplot);

NR = length(FWB);
jr1=1;
jr2=NR;
if ~isempty(nrplot) & nrplot>0
  jr1=nrplot;
  jr2=nrplot;
end


[mm,nn]=size(HH);

figure(nf);
clf;
contour(HH,[0 0],'k','Linewidth',1.6);
hold on;
contour(HH,[-5000:500:-10],'Color',[0.6 0.6 0.6]);
for jrg=jr1:jr2
  IJs=FWB(jrg).IJs;
  ns=size(IJs,1);
  for in=1:ns-1
    i1=IJs(in,1);
    j1=IJs(in,2);
    i2=IJs(in+1,1);
    j2=IJs(in+1,2);
    plot([i1 i2],[j1 j2],'r-','Linewidth',2);
  end
end

axis('equal');
set(gca,'xlim',[1 nn],...
	'ylim',[1 mm],...
	'xtick',[],...
	'ytick',[]);

nf=10;
dlmb=20;
dphi=10;
clr=[0.8 0.8 0.8];
plot_gridlines(dlmb,dphi, nf, clr, LON, LAT)
%bottom_text(txtb);



return