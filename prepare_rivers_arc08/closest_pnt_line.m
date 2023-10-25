    function [xm,ym,DM] = closest_pnt_line(A,B,C);
%
% Find a closest to C point on AB
%
% A(x,y) *--------------* B  (x,y)
%
%
%            C(x,y) *
%
% (1) I parameterize the line with t [0,1]
% x=(xb-xa)*t+xa=ax*t+xa
% y=(yb-ya)*t+ya=ay*t+ya
% (2) then I expressed distance between C and any pnt on AB
% in terms of t
% (3) Solved derivative dD/dt=0 -> found t
% (4) converted back to x and y
xa=A(1);
ya=A(2);
xb=B(1);
yb=B(2);
xc=C(1);
yc=C(2);

ax=xb-xa;
ay=yb-ya;
a=xc-xa;
b=yc-ya;
t=(a*ax+b*ay)/(ax^2+ay^2);

% if function does not have local min/max
% then t can be out of  bound
t=min([t,1]);
t=max([t,0]);

xm=ax*t+xa;
ym=ay*t+ya;

DA=sqrt((xa-xc).^2+(ya-yc).^2);
DB=sqrt((xb-xc).^2+(yb-yc).^2);
DM=sqrt((xm-xc).^2+(ym-yc).^2);

if DM>DA & DM>DB
  error('closest_pnt_line: Failed to find closest point');
end


%t=[0:0.01:1];
%D=sqrt((a-ax*t).^2+(b-ay*t).^2);




