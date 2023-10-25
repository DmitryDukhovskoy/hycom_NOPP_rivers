function sub_plot_xsectTSdecadal_hycom(ToPlot,Depth,...
	elonA,alatA,Fld,IJ,c1,c2,pthfig,id,sfig)
%
% ToPlot: variable whose cross-section this is
% Depth: depths
% Field name of variable with units
% yr1 and yr2: first and last year
% IJ number of points from left to right
% c1, c2: contour levels

fprintf('Plotting %s\n',Fld);

%nint=19;
%cnt=(c1:(c2-c1)/(nint):c2);
%cmp=colormap(jet(length(cnt)));

%cmp=[];
%nsint2=4;
%jjs=clrmp_zebra(cnt,cmp,nsint2);
%nsint2=1;

switch(Fld),
 case('T')
  nint= 240;
  CMP = create_colormap2_3(nint,c1,c2);
  cmp0= CMP.colormap;
  nav = 10;
  cmp = smooth_colormap(cmp0,nav);
  cnt = CMP.intervals;
 case('S')
  nint= 240;
%  CMP = create_colormap2_1(nint,c1,c2);
  CMP = create_colormap2_3(nint,c1,c2);
  cmp0= CMP.colormap;
  nav = 10;
  cmp = smooth_colormap(cmp0,nav);
  cnt = CMP.intervals;
end  




%xx=(1:length(IJ));
x1=elonA(IJ(1,2),IJ(1,1));
y1=alatA(IJ(1,2),IJ(1,1));
for ii=1:length(IJ)
  i=IJ(ii,1);
  j=IJ(ii,2);
  x2=elonA(j,i);
  y2=alatA(j,i);
  dx=distance_spheric_coord(y1,x1,y2,x2)*1e-3;  % spheric distance, km
  x1=x2;
  y1=y2;
  if (ii>1)
    xx(1:78,ii)=xx(:,ii-1)+dx;
  else
    xx(1:78,ii)=0;
  end
end;


nf=1;
figure(nf); clf;

pcolor(xx,Depth,ToPlot);
colormap(cmp)
shading interp;
caxis([c1 c2]);


%colormap(jjs);
xl1=0;
xl2=max(max(xx));

%set(gca,'ylim',[-200 0]);
set(gca,'ylim',[-4300 0],...
	'ytick',[-4500:500:0],...
	'xlim',[xl1 xl2],...
	'xtick',[0:500:xl2],...
	'tickdir','out',...
	'Color',[0 0 0]);

%keyboard


hght=[];
lngth=[];
mint=1;
mbx=mint;
fsz=11;
bxc='k';
posc=[0.9 0.06 0.9 0.06];
aend=0;
%  [az,axc]  = colorbar_horiz (cmp,cnt,hght,lngth,mint,fsz,bxc,posc,mbx,aend);
%[az,axc]  = colorbar_vert (cmp,cnt,hght,lngth,mint,fsz,bxc,posc,mbx,aend);
%colorbar_vert
colorbar

HT=title(['GDEM4, ',Fld]);
set(HT,'FontSize',11);

hold on;

txtb = 'sub_plot_xsectTSdecadal_gdem.m';
bottom_text(txtb,'pwd',1);

%contour(xx2,Depth,ToPlot,[-0.75 -0.75],'k','LineWidth',1);
%switch(Fld)
% case('T')
%  contour(xx,Depth,ToPlot,[0 0],'k','LineWidth',1);
% case('S')
%  contour(xx,Depth,ToPlot,[34.8 34.8],'k','LineWidth',1);
%end


if sfig==1
  fgnm=sprintf('%sGDEM4_xsect_%s_%s_shelf.png',pthfig,id,Fld);
  fprintf('Saving %s\n',fgnm);
  print('-dpng','-r250',fgnm)
end;
