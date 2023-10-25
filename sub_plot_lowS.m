function sub_plot_lowS(Fld,lTr,nf,HH,xlim1,xlim2,ylim1,ylim2,LON,LAT,stl);
% Plot low S on the shelf
hmsk=HH;
hmsk(HH<0)=nan;

switch(Fld)
 case('salin');
  c1=16.5;
  c2=35.5;
%  f_cmp=3;
  f_cmp=4;
%  cr1 = 30;
%  cr2 = 34;
%  dc  = 1;
end


% Colormap:
switch(f_cmp)
 case(1)
  na=50;
  cl1=colormap_gray(na);
  cl2=colormap_purple(na);
  cl3=colormap_blue(na);
  cl4=colormap_green(na);
  cl5=colormap_yellow(na);
  cl6=colormap_red(na);
  cmp=[cl2;cl3;cl4;cl5;cl6];
  cmp=smooth_colormap(cmp,10);
  nint=length(cmp);
  cnt=(c1:(c2-c1)/nint:c2);  % construct array of intervals
 case(2)
  CMP = colormap_mld1(200,c1,c2);
  cmp = CMP.colormap;
 case(3)
  CMP = colormap_sclr1(200,c1,c2);
  cmp = CMP.colormap;
 case(4)
  cmp = flipud(colormap_cold(200));
%  cmp = colormap_cold(200);
end

nint=length(cmp);


% Normalize:
%mT=max(max(Tr));
%Tr=Tr./mT*100;
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
%contour(lTr,[cr1:dc:cr2],'Color',[0.8 0.8 0.8]);
%contour(HH,[0 0],'k','Linewidth',1);
%contour(HH,[-5000:1000:-100],'Color',[0.7 0.7 0.7],'Linewidth',1);
contour(HH,[-20 -20],'Color',[0.7 0.7 0.7],'Linewidth',1);
contour(HH,[-100 -100],'Color',[0.7 0.7 0.7],'Linewidth',1);
contour(HH,[-500 -500],'Color',[0.7 0.7 0.7],'Linewidth',1);
contour(HH,[-2000 -2000],'Color',[0.7 0.7 0.7],'Linewidth',1);
contour(HH,[-3000 -3000],'Color',[0.7 0.7 0.7],'Linewidth',1);
%contour(HH,[-50:10:-5],'Color',[0.85 0.85 0.85],'Linewidth',1);
%caxis([0 2]);
caxis([c1 c2]);
colormap(cmp);
axis('equal');
set(gca,'xlim',[xlim1 xlim2],...
	'ylim',[ylim1 ylim2],...
        'xtick',[],...
	'ytick',[],...
	'Color',[0.9 0.9 0.9]);
clr=[0.95 0.95 0.95];
plot_gridlines(45,10,1,clr,LON,LAT);
title(stl,'Fontsize',12,'Interpreter','none');

cb = colorbar;
set(cb,'Fontsize',12,...
       'Position',[0.91 0.18 0.017 0.7],...
       'TickLength',0.018);
%keyboard


%freezeColors;
%pcolor(hmsk); shading flat;
%colormap([0 0 0]);




return