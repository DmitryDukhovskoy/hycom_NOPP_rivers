function sub_plotWnd_chck(sr,ur,vr,nn,mm,s,uu,vv,...
			  LN,LT,HH,INDX,DV,mc,nc);
% Plot interpoalted winds (sr,ur,vr) on ARCc0.72
% from CFSR/CFS fields(s)
% to check rotated vectors etc.
figure(10); clf;
pcolor(sr); shading flat;
caxis([-6 6]);
colorbar

hold on;
dv=7;
cf=0.5;
beta=20;
v_col=[.3 .3 .3];
lwd=1;
scl=2;
for iv=1:dv:nn
  for jv=1:dv:mm
    x1=iv;
    y1=jv;
    x2=x1+ur(jv,iv)*scl;
    y2=y1+vr(jv,iv)*scl;
     draw_arrowF(x1,x2,y1,y2,cf,beta,v_col,lwd);
  end
end
contour(LN,[-180:90:180],'color',[0.5 0.5 0.5]);
contour(LT,[60:10:89],'color',[0.5 0.5 0.5]);
contour(HH,[0 0],'c');
axis('equal');
stl=sprintf('%i/%2.2i/%2.2i',DV(1:3));
title(stl);

[a1,a2]=size(uu);
figure(12); clf;
pcolor(s); shading flat;
colorbar
hold on
caxis([-6 6]);

II=INDX.II;

for im=1:dv:nn
  for jm=1:dv:mm
%      for jk=1:dv:nI
  j0=II(jm,im);
  [jv,iv]=ind2sub([mc,nc],j0);
    x1=iv;
    y1=jv;
    x2=x1+uu(jv,iv)*scl;
    y2=y1+vv(jv,iv)*scl;
     draw_arrowF(x1,x2,y1,y2,cf,beta,v_col,lwd);
  end
end
axis('equal');
set(gca,'ylim',[410 a1]);


return