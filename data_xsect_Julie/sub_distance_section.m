% Calculate distance along section for small segments
% Distance normalize and scale to match total length of the section (Ltot)
% Xl, Yl - geogr coord of small segments along the section lines
%
function XL = sub_distance_section(Xl, Yl, Ltot);
XL0 = [];
XLd = [];
nn=length(Xl);
x1=Xl(1);
y1=Yl(1);
for ii=1:nn
  x2=Xl(ii);
  y2=Yl(ii);
% Distance from point 0:
  dxx=distance_spheric_coord(y1,x1,y2,x2)*1e-3;  % spheric distance, km
  XL0(ii) = distance_spheric_coord(Yl(1),Xl(1),y2,x2)*1e-3;  % spheric distance, from pnt0
%  XL(ii) = dxx;
  dxx=max([dxx,0.01]);
  dXX(ii)=dxx;
  x1=x2;
  y1=y2;
  XLd(ii)=sum(dXX(1:ii));
end
  
% Normalize and rescale to actual distance from point 0:
XL=XLd/max(XLd)*Ltot;


return

