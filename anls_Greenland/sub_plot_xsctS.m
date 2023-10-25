function sub_plot_xsctS(Z,X,S,U,Hb,stl,nfg,xl1,xl2,yl1,yl2);
% S vertical section 
%
figure(nfg);clf;

cs1=33.9;
cs2=34.9;
CMP = colormap_sclr2(400,cs1,cs2);
%  cmpS=CMP.colormap;
cmpS=colormap(parula(400));

figure(nfg); clf;
axes('Position',[0.09 0.4 0.81 0.5]);
pcolor(X,Z,S); shading interp;
colormap(cmpS);
caxis([cs1 cs2]);
hold on;
%contour(X,Z,S,[33.4:0.2:34.9],'k');
%contour(X,Z,S,[34 34],'k','linewidth',1.6);
contour(X,Z,U,[-1:0.1:0],'w');
contour(X,Z,U,[0 0],'w','linewidth',1.6);
%keyboard

hbx=[X(1);X;X(end)];
hby=[-5000;Hb;-5000];
fill(hbx,hby,[0 0 0]);


%xl1=0;
%xl2=max(X);
%yl1=-900;
%yl2=0;
set(gca,'tickdir','out',...
	'xlim',[xl1 xl2],...
	'ylim',[yl1 yl2],...
	'xtick',[0:20:xl2],...
	'ytick',[yl1:100:0]);
xlabel('Eastward Distance, km');

hb2=colorbar;
set(hb2,'Position',[0.91 0.4 0.012 0.5],...
 'Ticklength',0.025,...
 'Fontsize',12);
title(stl,'Interpreter','none');


return
