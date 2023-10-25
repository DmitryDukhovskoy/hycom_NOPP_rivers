function sub_plot_xsctT(Z,X,T,U,Hb,stl,nfg,xl1,xl2,yl1,yl2);
% T vertical section 
%
figure(nfg);clf;

ct1=0;
ct2=5;
CMP = create_colormap5(400,ct1,ct2);
cmpT=CMP.colormap;
%cmpS=colormap(parula(400));

figure(nfg); clf;
axes('Position',[0.09 0.4 0.81 0.5]);
pcolor(X,Z,T); shading interp;
colormap(cmpT);
caxis([ct1 ct2]);
hold on;
%%contour(X,Z,T,[0:0.5:5],'k');
%contour(X,Z,T,[-2:0.5:-0.1],'k--');
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
	'ytick',[yl1:100:0],...
        'Color',[0 0 0]);
xlabel('Eastward Distance, km');

hb2=colorbar;
set(hb2,'Position',[0.91 0.4 0.012 0.5],...
 'Ticklength',0.025,...
 'Fontsize',12);
title(stl,'Interpreter','none');


return
