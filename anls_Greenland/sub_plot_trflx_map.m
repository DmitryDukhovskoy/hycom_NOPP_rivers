function sub_plot_trflx_map(HH,LON,LAT,GC,hf1,...
			 fn,stl,xl1,xl2,yl1,yl2,c1,c2);

cl2 = colormap_red(100);
cl1 = colormap_blue(100);
for ik=1:2;
  cl2(ik,:) = [1 1 1];
  cl1(ik,:) = [1 1 1];
end
cl1 = flipud(cl1);
cmp = [cl1;cl2];
cmp = smooth_colormap(cmp,10);
ncnt = length(cmp);
%c1  = -2; % W/m
%c2  = 2;  
cntr= (c1:(c2-c1)/ncnt:c2);

Hg = HH;
Hg(1:380,:)=nan;
Hg(1100:end,:)=nan;
Hg(:,1050:end)=nan;
Hg(:,1:450)=nan;

fprintf('Plotting Greenland tracer flx contour ...\n');

domname = '0';
sub_plot_bath(Hg,LON,LAT,fn,domname);
%contour(Hg,[-100:10:-5],'Color',[0.85 0.85 0.85]);
%contour(Hg,[-450:50:-10],'Color',[0.6 0.6 0.6]);
%contour(Hg,[-4000:500:-50],'Color',[0.5 0.5 0.5]);
contour(Hg,[-5000:100:-50],'Color',[0.45 0.45 0.45]);
IIs = GC.cntr_Iindx;
JJs = GC.cntr_Jindx;
x   = GC.Distance_m*1e-3; % m->km

set(gca,'xtick',[],'ytick',[]);      

ngr = length(hf1);
for ik=1:ngr
  flx = hf1(ik);
  if isnan(flx); continue; end;
  icc  = max(find(cntr<=flx));
  if isempty(icc), icc=1; end;
  icc = min([icc, length(cmp)]);
  clr = cmp(icc,:);
  x0  = IIs(ik);
  y0  = JJs(ik);
%	scatter(x0,y0,8,clr,'o','filled');
  plot(x0,y0,'.','Markersize',26,'Color',clr);
end

% Plot 500km Markers
for km=0:500:max(x)
  d=abs(x-km);
  i0=find(d==min(d));
  if km==0
     plot(IIs(i0),JJs(i0),'k.','Markersize',15);
     plot(IIs(i0),JJs(i0),'kd','Markersize',7);
  else
    plot(IIs(i0),JJs(i0),'k.','Markersize',15);
  end
%  text(IIs(i0),JJs(i0),sprintf('%6.1f km',km),'Fontsize',14);
end

set(gca,'xlim',[xl1 xl2],...
	'ylim',[yl1 yl2],...
	'xtick',[],...
	'ytick',[]);

title(stl);

wdth=0.03;
hght=0.8;
mint=20;
fsz=12;
bxc=[0 0 0];
posc=[0.87 0.15 0.7 0.05]; 
mbx=mint;
aend=0;
[az,axc] = colorbar_vert(cmp,cntr,wdth,hght,mint,fsz,bxc,posc,mbx,aend);



return
