function [ic1,ic2,jc1,jc2] = sub_UVcntr(i0,j0,i1,j1,di,dj);
% Find segment of the grid point where U or V point is
% Note orientation of the segment
% depende on the general direction of the contour
%
%           * ic2, jc2
%           |
%  *(i1,j1) |       * (i0,j0)     Vertical segment contour is going up/down
%           |                     Flux is through vertical segment V*dy
%           * ic1,jc1             Here j1=j0 = jc
%
%
%                 * (i0,j0)
%    
%        *--------|--------* ic2, jc2 (can be ic1, jc1 - depends on contour direction)
% ic1,jc1         | V
%
%                  
%                 * (i1,j1)
%
ic=0.5*(i0+i1);
jc=0.5*(j0+j1);

if ic==i0
	ic1=ic-0.5;
	ic2=ic+0.5;
	jc1=jc;
	jc2=jc;
else
	ic1=ic;
	ic2=ic;
	jc1=jc-0.5;
	jc2=jc+0.5;
end

% !st point - no previous point
%if isempty(ic_old),
%  ic_old=ic1;
%  jc_old=jc1;
%end

if di<0 & ic2>ic1
  dmm=ic1;
  ic1=ic2;
  ic2=dmm;
elseif di>0 & ic1>ic2
  dmm=ic1;
  ic1=ic2;
  ic2=dmm;
end

if dj<0 & jc2>jc1
  dmm=jc1;
  jc1=jc2;
  jc2=dmm;
elseif dj>0 & jc1>jc2
  dmm=jc1;
  jc1=jc2;
  jc2=dmm;
end


return
