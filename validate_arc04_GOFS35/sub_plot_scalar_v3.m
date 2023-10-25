function sub_plot_scalar_v3(lTr,nf,HH,xlim1,xlim2,ylim1,ylim2,LON,LAT,stl,c1,c2,f_cmp);

hmsk=HH;
hmsk(HH<0)=nan;

% Colormap:
%c1 = 0;
%c2 = 0.2;
f_cmp =1;


nint = 200;
if f_cmp==1
%CMP = colormap_sclr1(nint,c1,c2);
CMP = create_colormap4(nint,c1,c2);
cmp = CMP.colormap;
cnt = CMP.intervals;
nint=length(cmp);
else
  cmp=colormap(parula(nint));
end


if (nf<0),
  nf=abs(nf);
  figure('Visible','off'); clf;
  fprintf('Figure window is turned off\n');
else  
  figure(nf);
  clf
end

pcolor(lTr); shading flat;
hold on;
%contour(HH,[0 0],'k','Linewidth',1);
contour(HH,[-5000:1000:-100],'Color',[0.7 0.7 0.7],'Linewidth',1);
%caxis([0 2]);
caxis([c1 c2]);
colormap(cmp);
axis('equal');
set(gca,'xlim',[xlim1 xlim2],'ylim',[ylim1 ylim2]);
set(gca,'xtick',[],'ytick',[]);
clr=[0.9 0.9 0.9];
plot_gridlines(45,10,1,clr,LON,LAT);
title(stl,'Fontsize',12,'Interpreter','none');

hght=[];
lngth=[];
mint=20;
mbx=mint;
fsz=12;
bxc='k';
posc=[0.81 0.1 0.8 0.06];
aend=0;
[az,axc]  = colorbar_vert(cmp,cnt,hght,lngth,mint,fsz,bxc,posc,mbx,aend);

freezeColors;
pcolor(hmsk); shading flat;
colormap([0 0 0]);




return