% Calculate distance along section for small segments 
% POP force to match HYCOM segments
% Distance normalize and scale to match total length of the section (Ltot)
% Xl, Yl - geogr coord of small segments along the section lines
%
function XL = sub_scale_POPdist(Xl, Yl, Chck);

npp = size(Chck,1);
XL  = 0;
XL0 = 0;
nn=length(Xl);

iS = 1;
for ik = 2:npp
  yy0 = Chck(ik,4);
  xx0 = Chck(ik,3);
  Lsgm = Chck(ik,5) - Chck(ik-1,5);

  dd = distance_spheric_coord(Yl,Xl,yy0,xx0)*1e-3;
  iE = find(dd==min(dd));

  x1=Xl(iS);
  y1=Yl(iS);
  icc = 0;
  XLd = [];
  for ii=iS+1:iE
    x2=Xl(ii);
    y2=Yl(ii);
  % Distance from point 0:
    dxx=distance_spheric_coord(y1,x1,y2,x2)*1e-3;  % spheric distance, km
    XL0(ii) = distance_spheric_coord(Yl(1),Xl(1),y2,x2)*1e-3;  % spheric distance, from pnt0
  %  XL(ii) = dxx;
    dxx=max([dxx,0.01]);
    icc = icc+1;
    dXX(icc)=dxx;
    x1=x2;
    y1=y2;
    XLd(icc)=sum(dXX(1:icc));
  end
  
%  keyboard  
% Normalize and rescale to actual distance from point 0:
  XLd = XLd/max(XLd)*Lsgm;
  XL(iS+1:iE) = XL(iS) + XLd;

  iS = iE;
end

XL = XL(:);

return

