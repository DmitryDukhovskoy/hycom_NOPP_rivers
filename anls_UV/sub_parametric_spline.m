function [Xi,Yi] = sub_parametric_spline(X,Y,nip);
% Interpolation of a 2D curve 
% represented as a paramteric curve
% X=x(t), Y=y(t)
% Here Parameter t is a cumulative distance or length of segments
% and the spline is the cumulative length spline
% See text book Quarteroni, p. 359
t(1)=0;
for i=1:length(X)-1
  L=sqrt((X(i+1)-X(i)).^2+(Y(i+1)-Y(i)).^2);
  t(i+1)=t(i)+L;
end
nt=length(t);
if isempty(nip), nip=4000; end;
%nip=4000;   % # of interpolation nodes
zip=[t(1):(t(nt)-t(1))/nip:t(nt)]; % interpolation points
Xi = spline(t,X,zip);
Yi = spline(t,Y,zip);

return