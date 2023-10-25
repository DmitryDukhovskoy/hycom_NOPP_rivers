function UVC=sub_UVpnts_contour(I,J,nin,HH);
% This is for NOT COLLOCATED U,V components !!!!
% Check if the components have been collocated
% in the archv files
% 
% Given a section/ contour I,J that goes through p-points
% Find corresponding contour going through U/V pnts
% used for flux calculation
%
% Collocate U,V,T,S variables
% along section or contour (I,J)
% in HYCOM grid
% nin - specify norm directed inside the contour (=1) or outside (=-1)
%     - if nin>1 - then this is an index that 
%                  indicates positive direction of the norm
%                  along the section
%                  
%  nin=1 or -1 is better to use for the contours/closed polygons
%              or group of connected sections that can form
%              a closed polygon by connecting the very 1st and the very last end points
%              nin=1 - means flows into this polygon is +
%
%  nin=index >1 - better to use for sections where inside/outside
%                cannot be properly defined
%         be cautious using index for polygons
%         as some segments on the contours may have positive direction
%         towards this point, when pointing outside the contour
%
%  
%      =[] - skip norm - not recommended for sections/contours
%            with "stepping" segments as along the steps
%            positive U or V component should in fact give
%            a negative flux !!!
%
% Return:  segments going through U/V pnts oriented from start pnt to end
%          p-point indices to interpolate dH and T, S for flux calculations (used in
%            sub_collocate2uv.m)
%          norm vectors at each segment (if nin - interior point with + flux is specified)
%
% The contour always goes through u(i,j) or v(i,j) points
% So that the direction you follow the contour doesn't matter
% The problem is "corner" points: in some cases (depending on the
% angle the contour changes) both u(i,j) and v(i,j) have to be
% added (to close off the contour), in others none are needed
% (interior corner grid cell)
% This is determined based on the orientation of P-vector (local
% direction) and u,v points wrt p-vector
%
%
%  HYCOM grid:
%
%          ------------   | V(i,j+1) ------------------
%          |                             |
%          |                             |
%        --- U(i,j)      *P(i,j)        ---  U(i+1,j)
%
%          |                             |
%          |                             |
%          ------------   | V(i,j) ------------------
%
%  Section follows U-V points
%  not centroids (P-points)
%   
%   
btx='sub_UVpnts_contour.m';

fprintf('%s\n',btx);


UVC = struct;
UVC.Norm=[];


ux=[1,0];
uy=[0,1];

ic1=0;
jc1=0;
ic2=0;
jc2=0;

nsc=length(I);
cc=0;
nch=round(0.2*nsc);
for is=1:nsc
  if mod(is,nch)==0, fprintf(' %5.2f done ...\n',is/nsc*100); end

  i0=I(is);
  j0=J(is);
  if is<nsc
    ip1=I(is+1);
    jp1=J(is+1);
    di=I(is+1)-I(is);
    dj=J(is+1)-J(is);
  else
    di=I(is)-I(is-1);
    dj=J(is)-J(is-1);
    ip1=I(is)+di;
    jp1=J(is)+dj;
  end
  if is>1
    im1=I(is-1);
    jm1=J(is-1);
  else
    im1=I(is)-di;
    jm1=J(is)-dj;
  end
%
% Checking correct segment line coodrinates: no points should be repeated
  if di==0 & dj==0
    fprintf('sub_UVpnts_contour: Line segm=%i repeated with previous point di/dj=0\n',is);
    keyboard
  end
%
  
% Get local direction pn vector and 
% un orientation wrt to pn
  pnx=ip1-i0;
  pny=jp1-j0;
  pn=[pnx;pny;0]./(sqrt(pnx^2+pny^2));

% un = u or v component:
% en is vector pointing towards U or V component
% from the p-point wrt to p-vector direction (contour)
  if pnx==0   % vertical segment
    en=[-1;0;0];
  else
    en=[0;-1;0];
  end
  
  en_pn=(en(1)*pn(2)-en(2)*pn(1));
  
% Same for the previous step:
  pnxL=i0-im1;
  pnyL=j0-jm1;
  if is==1
    en_pnL=en_pn;
    pnL=pn;
  else
%    pnxL=i0-im1;
%    pnyL=j0-jm1;
    pnL=[pnxL;pnyL;0]./(sqrt(pnxL^2+pnyL^2));
    if pnxL==0
      enL=[-1;0;0];
    else
      enL=[0;-1;0];
    end
  
    en_pnL=(enL(1)*pnL(2)-enL(2)*pnL(1));
  end
%
% Has the local orientation of the contour changed?
  pn_pnL=(pn(1)*pnL(2)-pn(2)*pnL(1));
  
%  keyboard
  
  
% Decide which components is needed:
% on a straight segments - only 1 components needed
%   section goes either through u or v-point
%   corner points - open corner - both components needed
%   closed corner - none
% u,v, both or none - depending
% on the direction of the contour & corner points
% when the contour turns, the uv-contour may cross
% the p-contour (going through p-points)
% this results in switching of the en-vector direction
% that point local direction from u/v point to the contour
% going through the p-points
  if en_pn==en_pnL & pn_pnL==0  % nothing changed, keep going
    if pny==0  % horizontal segment
      j1=j0-1; % grid pnt where V is coming from/to
      i1=i0;
    else
      j1=j0;
      i1=i0-1; % grid pnt where U is coming from/to
    end
    
    cc=cc+1;
    [ic1,ic2,jc1,jc2] = sub_UVcntr(i0,j0,i1,j1,di,dj);



    UVC.sgmX(cc,:)=[ic1,ic2];
    UVC.sgmY(cc,:)=[jc1,jc2];
    UVC.Ictr(cc)=is;
    UVC.gridIndx_IJ(cc,:)=[i0,j0];
    UVC.adjIndx_I1J1(cc,:)=[i1,j1];

%  ==============    
    
% Corner points:
% contour changed direction, U/V still on same side
% corner point - check both components to close the contour
%  Example of an "open" corner point:
%  |________________|
%  |                |
%  |                |
%  -       * p-point| 
%  |                |
%  |                |
%  |                |
%  |_______|________|___________
%  |       ^        |
%  | en    | pn-vect|   pn-vector indicates direction following the contour p-pnts
%  - <--   * p-point|  <--* p-pnt contour
%  |       |        |     |  
%  |       |        |     | en-vector pointing towards V-norm component
%  |       v en     |     v  on the u.v contour
%  |                |        here only 1 component needed to close contour
%  |                |
%  |_______|________|_____|______  u/v pnt contour
%         V comp          V comp
%
% Note direction of pn (clockwise or cntr/clckwise) is unimportant
% for detecting the segments and norms, however
% chosen direction probably should be consistent for 1 contour
%
  elseif en_pn==en_pnL & pn_pnL~=0
    % Check u-comp if needed
    % if u-comp at (i0,j0) points to a contour p-point - then not needed
    % closed corner
    i1=i0-1;
    j1=j0;
    du=min(sqrt((I-i1).^2+(J-j1).^2));
    % if v-comp at (i0,j0) points to a contour p-point - then not needed
    % closed corner
    i1=i0;
    j1=j0-1;
    dv=min(sqrt((I-i1).^2+(J-j1).^2));
% Sanity check: this cannot be, violation of a closed corner point
    if (du==0 & dv>0) | (du>0 & dv==0)
      fprintf('***ERR: collocate_UVC_sections: incosistency corner pnts\n');
      keyboard;
    end
    
    if du==0 & dv==0 % "closed" corner point, no U/V transport 
      continue;
    end
    
% For open corner points
% need to decide which side close first
% to make contour continuously connected
% with the previous segment
    if du>0 & dv>0  % both components are needed, "open" corner point
% U-normal
      j1=j0;
      i1=i0-1;
      
      ic1=0.5*(i1+i0);
      jc1=j0;
      if cc>1
        dsgm=min(sqrt((UVC.sgmX(cc,:)-ic1).^2+(UVC.sgmY(cc,:)-jc1).^2));
      else
	dsgm=0;
      end
      
      if dsgm<1
        cE=0;
      else
	cE=1;
      end
      cc=cc+1;
      
      [ic1,ic2,jc1,jc2] = sub_UVcntr(i0,j0,i1,j1,di,dj);
      UVC.sgmX(cc+cE,:)=[ic1,ic2];
      UVC.sgmY(cc+cE,:)=[jc1,jc2];
      UVC.Ictr(cc+cE)=is;
      UVC.gridIndx_IJ(cc,:)=[i0,j0];
      UVC.adjIndx_I1J1(cc,:)=[i1,j1];

% ===================
      
% V-normal      
      j1=j0-1;
      i1=i0;
      
      cE=-cE; % adjust cc slot if needed, either 0 or -1
      
      cc=cc+1;

      [ic1,ic2,jc1,jc2] = sub_UVcntr(i0,j0,i1,j1,di,dj);
      UVC.sgmX(cc+cE,:)=[ic1,ic2];
      UVC.sgmY(cc+cE,:)=[jc1,jc2];
      UVC.Ictr(cc+cE)=is;
      UVC.gridIndx_IJ(cc,:)=[i0,j0];
      UVC.adjIndx_I1J1(cc,:)=[i1,j1];
      
% ===================
    end % if du>0
    
      
% Complicated case when the contour crosses
% the line where U/V transport is actually estimated
% This is reflected in change of orientation of en-vector wrt pn
% i.e. normal component from left-side switches to right-sight
% wrt to the local contour direction
% In this case choose normal component that 
% points along the contour (either backward or fwd)
  elseif en_pn~=en_pnL & pn_pnL~=0
    i1=i0-1;
    j1=j0;
    du=min(sqrt((I-i1).^2+(J-j1).^2));
    i1=i0;
    j1=j0-1;
    dv=min(sqrt((I-i1).^2+(J-j1).^2));
% Sanity check: this cannot be:    
    if (du==0 & dv==0) | (du>0 & dv>0)
      fprintf('***ERR: collocate_UVC_sections: 2. Incosistency corner pnts\n');
      keyboard;
    end

    if du==0 % need u-normal 
      j1=j0;
      i1=i0-1;
    else   % v-normal
      j1=j0-1;
      i1=i0;
    end

    cc=cc+1;

    [ic1,ic2,jc1,jc2] = sub_UVcntr(i0,j0,i1,j1,di,dj);
    UVC.sgmX(cc,:)=[ic1,ic2];
    UVC.sgmY(cc,:)=[jc1,jc2];
    UVC.Ictr(cc)=is;
    UVC.gridIndx_IJ(cc,:)=[i0,j0];
    UVC.adjIndx_I1J1(cc,:)=[i1,j1];
% ===================

  end  % if
  
end  % sect

% 
%
% Find norms if needed (for simple sections = straight lines  -------
% of vertical lines
% not needed and +/- X-Y orientation is taken) 
% Norms are needed for sections with "zig-zags" when
% in/out flow is 
% this is the case
% for polygons/contours/ boxes
% Norms are
% directed inside or outisde the contour
% First, construct a contour where fluxes 
% are calculated
% Second, go around the contour and 
% determine the direction of the norm 
% the end of the normal vector should be inside /outisde
% of the contour wrt to the centroid
% 

lpnt = 0;
if nin>1
  [jP,iP]=ind2sub(size(HH),nin);
  lpnt=1;
end

  
sgmx=UVC.sgmX;
sgmy=UVC.sgmY;

if ~isempty(nin); 

  % Combine contour going through u/v points
  ns=length(sgmx);
  Iuv=[];
  Juv=[];
  if sgmx(1,1)<sgmx(end,1)
    x1=min(sgmx(1,:));
    x2=max(sgmx(end,:));
  else
    x1=max(sgmx(1,:));
    x2=min(sgmx(end,:));
  end

  if sgmy(1,1)<sgmy(end,1)
    y1=min(sgmy(1,:));
    y2=max(sgmy(end,:));
  else
    y1=max(sgmy(1,:));
    y2=min(sgmy(end,:));
  end

  Iuv(1,1)=x1;
  Juv(1,1)=y1;
  xE=x2;
  yE=y2;
  for ik=1:ns
    x1=sgmx(ik,1);
    x2=sgmx(ik,2);
    y1=sgmy(ik,1);
    y2=sgmy(ik,2);
    Iuv=[Iuv;0.5*(x1+x2)];
    Juv=[Juv;0.5*(y1+y2)];
  end
  Iuv(end+1)=xE;
  Juv(end+1)=yE;
%  x1=xE+2*(sgmx(1,end)-sgmx(1,end-1)); % extend segment to guarantee norm is in at the endpoints
%  y1=sgmy(end,1)+2*(sgmy(1,end)-sgmy(1,end-1));
%  Iuv(ns,1)=x1;
%  Juv(ns,1)=y1;
  if nin>1, % for sections, create a polygon, inside - norm is >0
    Iuv(end+1)=iP;
    Juv(end+1)=jP;
  end
  Iuv(end+1)=Iuv(1);
  Juv(end+1)=Juv(1);
  
%  keyboard

%CC = Centroid([I,J]);
  for ik=1:ns
    x1=sgmx(ik,1);
    x2=sgmx(ik,2);
    y1=sgmy(ik,1);
    y2=sgmy(ik,2);
    xs0=0.5*(x1+x2); % midpoint
    ys0=0.5*(y1+y2);
  % Rotate by 90 degrees:
  % around mid point
    aa=x2-xs0;
    bb=y2-ys0;
    rx=-bb;
    ry=aa;
    ll=sqrt(rx*rx+ry*ry);
    nx=rx./ll;
    ny=ry./ll;
    rxx=xs0+0.01*nx; % shorten norm to make it stay within the grid cell
    ryy=ys0+0.01*ny;

    [inn,onn]=inpolygon(rxx,ryy,Iuv,Juv);
    if onn, inn=~inn; end % point on the edge, happens at corner points
    if abs(nin==1)
      if (inn & nin==1) | ...
       	 (~inn & nin==-1)
	       UVC.Norm(ik,1)=nx;
	       UVC.Norm(ik,2)=ny;
      else
        UVC.Norm(ik,1)=-nx;
	       UVC.Norm(ik,2)=-ny;
      end
    else
      if inn 
       	UVC.Norm(ik,1)=nx;
	       UVC.Norm(ik,2)=ny;
      else
	       UVC.Norm(ik,1)=-nx;
       	UVC.Norm(ik,2)=-ny;
      end
    end
    
    
  %  keyboard
  end
  
  


end % nin

%
% Order segments so that their endpoints are connected
sgmx=UVC.sgmX;
sgmy=UVC.sgmY;

% find local direction at the start of the contour
% and orient the 1st segment
Is=UVC.gridIndx_IJ(:,1);
Js=UVC.gridIndx_IJ(:,2);
di1=Is(2)-Is(1);
dj1=Js(2)-Js(1);

ii=1;
i1=sgmx(ii,1);
i2=sgmx(ii,2);
j1=sgmy(ii,1);
j2=sgmy(ii,2);
if dj1==0 
  if di1<0
    io1=max([i1,i2]);
    io2=min([i1,i2]);
  else
    io1=min([i1,i2]);
    io2=max([i1,i2]);
  end
  sgmx(ii,1)=io1;
  sgmx(ii,2)=io2;
else
  if dj1<0
    jo1=max([j1,j2]);
    jo2=min([j1,j2]);
  else
    jo1=min([j1,j2]);
    jo2=max([j1,j2]);
  end
  sgmy(ii,1)=jo1;
  sgmy(ii,2)=jo2;
end

for ii=2:ns
  iold=sgmx(ii-1,2);
  jold=sgmy(ii-1,2);
		i1=sgmx(ii,1);
		i2=sgmx(ii,2);
		j1=sgmy(ii,1);
		j2=sgmy(ii,2);

  if (i1~=iold & i2~=iold) | ...
     (j1~=jold & j2~=jold)
    fprintf('sub_UV_pnts_contour: segments are not connected, ii=%i ...\n',ii);
    keyboard
  end

  if i1~=iold
    sgmx(ii,1)=i2;
    sgmx(ii,2)=i1;
  end
  if j1~=jold
    sgmy(ii,1)=j2;
    sgmy(ii,2)=j1;
  end
end

% Checking:
%  figure(11); clf; hold on;
iold=sgmx(1,1);
jold=sgmy(1,1);
for ii=1:ns
		i1=sgmx(ii,1);
		i2=sgmx(ii,2);
		j1=sgmy(ii,1);
		j2=sgmy(ii,2);
		
		di(ii)=sqrt((i1-iold).^2+(j1-jold).^2);
		iold=i2;
		jold=j2;
%  plot([iold i1],[jold j1],'r-');
end

if max(di)>0,
		fprintf('sub_UV_pnts_: segments failed to be correctly oriented to have end-point connection\n');
		keyboard;
end


UVC.sgmX=sgmx;
UVC.sgmY=sgmy;

% --------------------------------

f_chck=0;
if f_chck==1
  fprintf('Plotting: Checking contour and u/v norms\n');
  figure(10); clf;
  hold on;
  contour(HH,[0 0],'k');
  contour(HH,[-5000:1000:0],'Color',[0.6 0.6 0.6]);
  axis('equal');
%  set(gca,'xlim',[800 1150],...
%	  'ylim',[800 1100]);
  set(gca,'xlim',[min(sgmx(:,1))-40 max(sgmx(:,1))+40],...
	  'ylim',[min(sgmy(:,1))-40 max(sgmy(:,1))+40]);
  plot(I,J,'.-');
  for ik=1:cc
    x=UVC.sgmX(ik,:);
    y=UVC.sgmY(ik,:);
    plot(x,y,'r-');

%
% Plot norms if defined
    if ~isempty(nin)
      nx=UVC.Norm(ik,1);
      ny=UVC.Norm(ik,2);

      xs0=mean(x); % midpoint
      ys0=mean(y);

      rxx=xs0+nx;
      ryy=ys0+ny;
      plot([xs0 rxx],[ys0 ryy],'Color',[0 1 0.5]);
      
    end
  
  end
  
  bottom_text(btx,'pwd',1);
  
  keyboard
  
  
end

% CHeck segments


  
return
