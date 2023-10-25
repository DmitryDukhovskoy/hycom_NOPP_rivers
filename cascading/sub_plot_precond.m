function sub_plot_precond(A)

nf=A.nf;
F=A.fld;
fldnm=A.fldnm;
xlim1=A.xl1;
xlim2=A.xl2;
ylim1=A.yl1;
ylim2=A.yl2;
c1 = A.c1c2(1);
c2 = A.c1c2(2);
HH=A.HH;
LON=A.LON;
LAT=A.LAT;
stl=A.title;

fprintf('Plotting %s...\n',stl);
fprintf('%s\n',fldnm);

nint=200;
switch(fldnm)
 case('prob');
  CMP = colormap_YOR(nint,c1,c2);
  cmp = CMP.colormap;
  for k=1:5
    cmp(k,:)=[1 1 1];
  end
  cmp=smooth_colormap(cmp,15);
  cmp(1,:)=[1 1 1];
  cnt = CMP.intervals;
  nint=length(cmp);
 case('temp');
  CMP = colormap_PBYR(nint,c1,c2);
  cmp = CMP.colormap;
  cnt = CMP.intervals;
  nint=length(cmp);
 case('salin');
  CMP = colormap_sclr2(nint,c1,c2);
  cmp = CMP.colormap;
  cnt = CMP.intervals;
  nint=length(cmp);
  
end

hmsk=HH;
hmsk(HH<0)=nan;

figure(nf); clf;
pcolor(F); shading flat;
hold on;
%contour(HH,[0 0],'k','Linewidth',1);
contour(HH,[-5000:1000:-100],'Color',[0.7 0.7 0.7],'Linewidth',1);
contour(HH,[-500 -500],'Color',[0.5 0.5 0.5],'Linewidth',1.6);
%caxis([0 2]);
caxis([c1 c2]);
colormap(cmp);
axis('equal');
set(gca,'xlim',[xlim1 xlim2],'ylim',[ylim1 ylim2]);
set(gca,'xtick',[],'ytick',[]);

freezeColors;
pcolor(hmsk); shading flat;
colormap([0 0 0]);
freezeColors;

clr=[0.9 0.9 0.9];
plot_gridlines(45,10,1,clr,LON,LAT);

colormap(cmp);

title(stl,'Fontsize',14,'Interpreter','none');

hb=colorbar;
set(hb,'Position',[0.89 0.2 0.02 0.6],...
       'Fontsize',16,...
       'TickLength',0.022);

return