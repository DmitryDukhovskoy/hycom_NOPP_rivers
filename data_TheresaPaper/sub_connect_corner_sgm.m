Not finished does not correctly identifies corner points for all possible cases


function [Isct,Igrn] = sub_connect_corner_sgm(IIs,JJs,Ig,Jg,icst,jcst,Intgr,Intgt,xP,yP,isct,f_map);
%
% Input: IIs, JJs - gate section indices
%        Ig,Jg - whole Gr contour
%        icst, jcst - coast coordinate
%        Intgr - intercept index for Gr contour
%        Intgt - intercept index for gate section
%        xP,yP - check point in the box
%        isct - section 1 or 2 wrt to the contour direction, 
%               i.e. it means contour indices increase from 1 to 2
%
% Output: Isct - intercept/cornerindex of the gate section to close the box
%         Igrn  - interecpt index of the Greenland contour to close the box
%
% At the corner points: need to discard intercept point or include depend
% ing on the position of the contour and gate section wrt to each other and
% Greenlnd coast
%
%       |
%       |     Gr Shelf ->
%       |
%   U  --            *  (i,j)  - intersection point   <--------  Gate section going left from coast
%       |
%       |  
%       |              V(i,j)
%    --------------- | ------------
%       |
%       |                        
%       |
%   U  --            *  (i,j-1)  -  Box  
%       |
%       |            ^  Gr contour going up/down
%       |            |
%       |            |
%       |            |
%
%  Here: need to include (i,j) for V flux along the gate section and (i,j-1) for U flux along
%        the Grenlan contour to close the corner
%
% Approach: find where U,V contours following gate section & Greenl contour intercept and close
%    off the box in the corner
% 

% Determine the corner point wrt to the center of the box:
% NW - upper left, NE - upper right, SW - lower lest, SE - lower right
%

% Gr contour, local direction at the intercept
% take 1 point in the box and 1 at the corner
if isct==1 % Gr contour entering the box
		ig1=Ig(Intgr);  % 1 point before intercept, back in the contour
		jg1=Jg(Intgr);
		ig2=Ig(Intgr+1);    % interecpt
		jg2=Jg(Intgr+1);
else
		ig1=Ig(Intgr-1);  % 1 point before intercept, back in the contour
		jg1=Jg(Intgr-1);
		ig2=Ig(Intgr);    % interecpt
		jg2=Jg(Intgr);
end

di=ig2-ig1;
dj=jg2-jg1;
btt=atan2(dj,di);

if dj==0, % V flux
% start/end points of V segment at point ig1, hg1
  x1gbck=ig1-0.5*cos(btt);
  x1gfwd=ig1+0.5*cos(btt);
  y1gbck=jg1-0.5;
  y1gfwd=jg1-0.5;
% start/end points of V segm pnt ig2,jg2 - intercept
  x2gbck=ig2-0.5*cos(btt); % this should match x1gfwd
  x2gfwd=ig2+0.5*cos(btt); % 
  y2gbck=jg2-0.5;
  y2gfwd=jg2-0.5;
else   % U flux
% end/start pnt of U segment 
  x1gbck=ig1-0.5;
  x1gfwd=ig1-0.5;
  y1gbck=jg1-0.5*sin(btt);
  y1gfwd=jg1+0.5*sin(btt);
%
  x2gbck=ig2-0.5;
  x2gfwd=ig2-0.5;
  y2gbck=jg2-0.5*sin(btt);
  y2gfwd=jg2+0.5*sin(btt);
end

% Section end closest to the coast:
dd=sqrt((IIs-icst).^2+(JJs-jcst).^2);
if dd(1)<dd(end)
		st1=1; % Greenl coast at index 1 of the segment
else
		st1=0;
end

if st1==1   % coast at the beginning of the section
  i1=IIs(Intgt-1);
  j1=JJs(Intgt-1);
  i2=IIs(Intgt);
  j2=JJs(Intgt);

  di=i2-i1;
  dj=j2-j1;

  alf=atan2(dj,di);
  if dj==0  % V flux
    x1bck=i1-0.5*cos(alf);
    x1fwd=i1+0.5*cos(alf);
    y1bck=j1-0.5;
    y1fwd=j1-0.5;
    x2bck=i2-0.5*cos(alf);
    x2fwd=i2+0.5*cos(alf);
    y2bck=j2-0.5;
    y2fwd=j2-0.5;
  else   % Uflux
    x1bck=i1-0.5;
    x1fwd=i1-0.5;
    y1bck=j1-0.5*sin(alf);
    y1fwd=j1+0.5*sin(alf);
    x2bck=i2-0.5;
    x2fwd=i2-0.5;
    y2bck=j2-0.5*sin(alf);
    y2fwd=j2+0.5*sin(alf);
  end

else % coast at the end of the gate section
  i1=IIs(Intgt+1);
  j1=JJs(Intgt+1);
  i2=IIs(Intgt);
  j2=JJs(Intgt);

  di=i2-i1;
  dj=j2-j1;

  alf=atan2(dj,di);
  if dj==0  % V flux
    x1bck=i1-0.5*cos(alf);
    x1fwd=i1+0.5*cos(alf);
    y1bck=j1-0.5;
    y1fwd=j1-0.5;
    x2bck=i2-0.5*cos(alf);
    x2fwd=i2+0.5*cos(alf);
    y2bck=j2-0.5;
    y2fwd=j2-0.5;
  else   % Uflux
    x1bck=i1-0.5;
    x1fwd=i1-0.5;
    y1bck=j1-0.5*sin(alf);
    y1fwd=j1+0.5*sin(alf);
    x2bck=i2-0.5;
    x2bck=i2-0.5;
    y2bck=j2-0.5*sin(alf);
    y2fwd=j2+0.5*sin(alf);
  end

end  

%keyboard

%
% Select UV segments that close off the corner
% wrt to the interior point in the box
SGt=[x1bck, y1bck; ...
     x2bck, y2bck; ...
     x1fwd, y1fwd; ...
     x2fwd, y2fwd];

SGr=[x1gbck, y1gbck; ...
     x2gbck, y2gbck; ...
     x1gfwd, y1gfwd; ...
     x2gfwd, y2gfwd];
%
% 1) select segments that have intercept point Greenl cntr - gate section
%
%                      Gr contour going through u pnts
%                     | segm2
%                     |
%                     |    
%                <------  U       *(i,j)
%                     |
%  x1                 |               ^
%  y1 segm 1          |x2,y2          | 
%   ----------------------------------|------  Section going through V pnts
%                     |  segm 2 (section)   
%    Box pnt          |               
%     *               |               
%                     | segm 1 (Gr contour)
%                     |               
%                     |               
%
%  segm 1 (section) and segm 1(Gr contour) close the box and box pnt is inside
%
%

icls=0;
II=[];
for kk=1:4
  xgt=SGt(kk,1);
  ygt=SGt(kk,2);
  dd=sqrt((SGr(:,1)-xgt).^2+(SGr(:,2)-ygt).^2);
  J=find(dd==0);
  for ll=1:length(J)
% Check orientation wrt to gate contour and compare with box point
				igt=kk;
				if mod(igt,2)==0, igt=igt-1; end
				x1gt=SGt(igt,1);
				y1gt=SGt(igt,2);
				x2gt=SGt(igt+1,1);
				y2gt=SGt(igt+1,2);

    igr=J(ll);
				if mod(igr,2)==0, igr=igr-1; end
				x1gr=SGr(igr,1);
				y1gr=SGr(igr,2);
				x2gr=SGr(igr+1,1);
				y2gr=SGr(igr+1,2);

% Pick segment that is in the box wrt to Pnt
    D = sign(orientation(x1gt,y1gt,x2gt,y2gt,xP,yP));
    R1 = orientation(x1gt,y1gt,x2gt,y2gt,x1gr,y1gr);
    R2 = orientation(x1gt,y1gt,x2gt,y2gt,x2gr,y2gr);
%
% Correct segment should have 0 - (i.e being on the section segmn)
% and be on the same side of section segm as the box pnt    
    if (R1==0 | R2==0) & (R1==D | R2==D)
%
% Now check that gate segment is on the same side as pnt
% wrt to Gr contour segments
      D2 = sign(orientation(x1gr,y1gr,x2gr,y2gr,xP,yP));
      P1 = orientation(x1gr,y1gr,x2gr,y2gr,x1gt,y1gt);
      P2 = orientation(x1gr,y1gr,x2gr,y2gr,x2gt,y2gt);

      if P1==D2 | P2==D2  % should be only 1 case
        icls=icls+1;
        II(icls,1)=kk;
        II(icls,2)=J(ll);
      end
    end
  end
end



if icls==0
  fprintf('sub_connect_corner_sgm.m:\n');
  fprintf('Connecting corners: could not find UV segments to close the box\n');
  keyboard
end
if icls>1
  fprintf('sub_connect_corner_sgm.m:\n');
  fprintf('Connecting corners: Found >1  UV segment pairs to close the box\n');
  keyboard
end

igt=II(1);
if mod(igt,2)==0, igt=igt-1; end
x1gt=SGt(igt,1);
y1gt=SGt(igt,2);
x2gt=SGt(igt+1,1);
y2gt=SGt(igt+1,2);
% Retrieve corresponding p-point:
GtPnt=[min([x1gt,x2gt])+0.5,min([y1gt,y2gt])+0.5];
dd=sqrt((IIs-GtPnt(1)).^2+(JJs-GtPnt(2)).^2);
Isct=find(dd==min(dd));


igr=II(2);
if mod(igr,2)==0, igr=igr-1; end
x1gr=SGr(igr,1);
y1gr=SGr(igr,2);
x2gr=SGr(igr+1,1);
y2gr=SGr(igr+1,2);
GrPnt=[min([x1gr,x2gr])+0.5,min([y1gr,y2gr])+0.5];
dd=sqrt((Ig-GrPnt(1)).^2+(Jg-GrPnt(2)).^2);
Igrn=find(dd==min(dd));


% plot corner UV segments
% if debug mode
f_chck=1;
if f_chck==1 & f_map>0
%  figure(f_map); 
  plot([x1bck x2bck],[y1bck y2bck],'-','Color',[0.4 0.4 0.4]);
  plot([x1fwd x2fwd],[y1fwd y2fwd],'-','Color',[0.4 0.4 0.4]);
  plot([x1gbck x2gbck],[y1gbck y2gbck],'-','Color',[0.4 0.4 0.4]);
  plot([x1gfwd x2gfwd],[y1gfwd y2gfwd],'-','Color',[0.4 0.4 0.4]);
% Selected uv segments:
  plot([x1gr x2gr],[y1gr y2gr],'-','Color',[1 0.6 0]);
  plot(GrPnt(1),GrPnt(2),'.','Markersize',14,'Color',[1 0.6 0]);
  plot([x1gt x2gt],[y1gt y2gt],'-','Color',[0 0.9 0.7]);
  plot(GtPnt(1),GtPnt(2),'.','Markersize',14,'Color',[0 0.9 0.7]);
end


return
