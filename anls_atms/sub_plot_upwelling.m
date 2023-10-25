function sub_plot_upwelling(HH,LN,LT,fn,Mnrm,GRC,stl);
% Plot upwelling index
%
%DV = datevec(dnmb);
%yr = DV(1);
%im = DV(2);
ngr = length(GRC.X);


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
%c1  = -0.3; % m3/s
%c2  = 0.3;  
c1  = -0.2; % m3/s
c2  = 0.2;  
cntr= (c1:(c2-c1)/ncnt:c2);


xyl = [465, 355;...
       965, 1100];

sub_plot_bath2(HH,LN,LT,fn,xyl);
%      contour(HH,[-100:10:-5],'Color',[0.3 0.3 0.3]);
%      contour(HH,[-3500:100:-10],'Color',[0.6 0.6 0.6]);
      contour(HH,[-4000:500:-50],'Color',[0.8 0.8 0.8]);
set(gca,'xtick',[],'ytick',[]);      

for ik=1:ngr
  mnrm = Mnrm(ik);
  icc  = max(find(cntr<=mnrm));
  if isempty(icc), icc=1; end;
  icc = min([icc, length(cmp)]);
  clr = cmp(icc,:);
  x0  = GRC.Iocn(ik);
  y0  = GRC.Jocn(ik);
%	scatter(x0,y0,8,clr,'o','filled');
  plot(x0,y0,'.','Markersize',20,'Color',clr);
end

%stl = sprintf('Upwelling index, CFSR/CFSv2, %i/%2.2i, m3/s per 1m',...
%	      yr,im);
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