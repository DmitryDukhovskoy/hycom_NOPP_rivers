function sub_plot_hice(A,Lmsk,xl1,xl2,yl1,yl2,c1,c2,fn,cmp);
%
%
%Lmsk = A*0+1;
%Lmsk(isnan(A))=0;

%c1=0;
%c2=5;
if isempty(cmp)
nint=200;
CMP = colormap_WB(nint,c1,c2);
cmp=CMP.colormap;
cmp(1,:)=[1 1 1];
else 
  nint=length(cmp);
end;

%cmp=smooth_colormap(cmp,15);

lcmp=[0.5 0.5 0.5; 0.9 0.9 0.9];

figure(fn); clf;
set(gcf,'Position',[1209 280 1019 1055]);
axes('Position',[0.05 0.1 0.8 0.8]);
hold on;
pcolor(Lmsk); shading flat;
colormap(lcmp);
freezeColors;

pcolor(A); shading flat;
caxis([c1 c2]);
colormap(cmp);

axis('equal');

set(gca,'xlim',[xl1 xl2],...
        'ylim',[yl1 yl2],...
        'xtick',[],...
        'ytick',[]);
%        'Color',[0 0 0]);

hb=colorbar;
set(hb,'Position',[0.87 0.11 0.02 0.8],...
       'Ticklength',0.03,...
       'Fontsize',14);


return
