function sub_plot_hice(A,xl1,xl2,yl1,yl2);
%
%
Lmsk = A*0+1;
Lmsk(isnan(A))=0;

c1=0;
c2=5;
nint=200;
CMP = colormap_WB(nint,c1,c2);
cmp=CMP.colormap;
cmp(1,:)=[1 1 1];

figure(5); clf;
hold on;
pcolor(A); shading flat;
caxis([c1 c2]);
colormap(cmp);

axis('equal');

set(gca,'xlim',[xl1 xl2],...
        'ylim',[yl1 yl2],...
        'xtick',[],...
        'ytick',[],...
        'Color',[0 0 0]);

hb=colorbar;
set(hb,'Position',[0.82 0.11 0.02 0.8],...
       'Ticklength',0.03,...
       'Fontsize',12);





return
