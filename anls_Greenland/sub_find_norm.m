function [nx,ny]=sub_find_norm(II,JJ,IJC); 
% Find norm orientation for 
% segment with coord. II,JJ
% in the direction to the coast line IJC

x0=mean(II);
y0=mean(JJ);
X=IJC(:,1);
Y=IJC(:,2);

dd=sqrt((X-x0).^2+(Y-y0).^2);
i0=find(dd==min(dd),1);
x0c=X(i0);
y0c=Y(i0);

uu=x0c-x0;
vv=y0c-y0;
U=sqrt(uu.^2+vv.^2);
%plot(x0,y0,'k*')
%plot(x0c,y0c,'g*')
nx=uu/U;
ny=vv/U;

%keyboard

return