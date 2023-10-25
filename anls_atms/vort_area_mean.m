    function [VORT,INDX]=vort_area_mean(X,Y,U,V,RR,INDX);
%
% Calculate area-mean vorticity
% using circulation theorem
% see Bourassa and Ford, 2010
%
% Input:
%   U,V - wind field
%   X,Y - geogr. coordinates of U,V
%         Note that U,V, X, Y has to be
%         "glued" along the 0 meridian
%         so that X(1)<X(2)... <X(N), where
%         X(1) is the westernmost location
%   RR - ring size ("diameter"), km, will be rounded 
%        to the nearest grid points
%
%   INDX - struct array of indices of ring contours
%
% Output: VORT - structured array
%
clear VRT AREA
l=min(size(X));
if l==1
  [X,Y]=meshgrid(X,Y);
end;

[m,n]=size(U);
[IIND,JIND]=meshgrid([1:n],[1:m]);

%
get_indx=logical(0);
if isempty(INDX)
  get_indx=logical(1);
  fprintf('Getting ring indices ...\n');
end

VRT=zeros(m,n);
AREA=zeros(m,n);

for i=1:n
%  fprintf('i=%i \n',i);
  for j=1:m
    lI=sub2ind([m,n],j,i);

    if isnan(X(j,i)), continue; end; % point to skip
    
    if get_indx
      INDX(lI).IJ=[];
      INDX(lI).XY=[];

      j0=j;
      i0=i;
      
      if j==m | i==n, continue; end;
      if j==1 | i==1, continue; end;
      
      ln0=X(j,i);
      lt0=Y(j,i);
      ln1i=X(j,i+1);
      lt1i=Y(j,i+1);
      ln1j=X(j+1,i);
      lt1j=Y(j+1,i);

      if isnan(ln0) | isnan(ln1i) | isnan(ln1j),
	continue;
      end
      
      dx=distance_spheric_coord(lt0,ln0,lt1i,ln1i)*1e-3;
      dy=distance_spheric_coord(lt0,ln0,lt1j,ln1j)*1e-3;
      if dx==0 | dy==0,
	error('Check dx, dy ...');
      end
      

      di=round(0.5*RR/dx);
      dj=round(0.5*RR/dy);
      if di==0 & RR/dx>0
        di=1;  % ring size is smaller than grid spcing
      end
      if dj==0 & RR/dy>0
        dj=1;  % ring size is smaller than grid spcing
      end
      

  % if ring size is smaller than horizontal spacing
  % stick with the spacing
  % Note, if di or dj =0 for given RR
  % the "ring" shape may not be regular
  % as dx or dy may be very different
  % from RR:
      if di==0, 
				di=max([1,di]);
				r=di*dx;
				dj=round(r/dy);
      elseif dj==0
				dj=max([1,dj]);
				r=dj*dy;
				di=round(r/dx);
      end

  % Construct indices:
      cc=1;
      sgn=complex(0,1);
      iv=complex(0,1);
      clear IJ
      diN=i+di;
      diN=min([diN,n]);
      IJ(1,1)=diN;
      IJ(1,2)=j;

      cch=0;
      while cch<1000
	cch=cch+1;
	iold=IJ(cc,1);
	jold=IJ(cc,2);
	inew=iold+real(sgn);
	jnew=jold+imag(sgn);
	dI=abs(inew-i);
	dJ=abs(jnew-j);
	if dI>di | dJ>dj | inew>n | inew<1 | jnew>m | jnew<1
	  sgn=sgn*iv;
	else
	  cc=cc+1;
	  IJ(cc,1)=inew;
	  IJ(cc,2)=jnew;
	end
	d0=sqrt((inew-IJ(1,1))^2+(jnew-IJ(1,2))^2);
	if d0==0, break; end;
  %      fprintf('sgn=%i %ii,inew=%i, jnew=%i, dI=%i, dJ=%i, cc=%i\n',...
  %	      real(sgn),imag(sgn),inew,jnew,dI,dJ,cc);
  %      keyboard
      end

      in=inpolygon(i,j,IJ(:,1),IJ(:,2));

      if d0~=0;
	fprintf('ERROR: i=%i j=%i\n',i,j);
	keyboard
	error('vort_area_mean: Could not locate indices ...');
      end

      if in~=1,
	error('vort_area_mean: Wrong indces - point is not inside ...');
      end

  % eliminate "Corners":
      d=sqrt((IJ(:,1)-(i-di)).^2+(IJ(:,2)-(j-dj)).^2);
      iM=find(d==0);
      IJ(iM,:)=nan;
      d=sqrt((IJ(:,1)-(i-di)).^2+(IJ(:,2)-(j+dj)).^2);
      iM=find(d==0);
      IJ(iM,:)=nan;
      d=sqrt((IJ(:,1)-(i+di)).^2+(IJ(:,2)-(j+dj)).^2);
      iM=find(d==0);
      IJ(iM,:)=nan;
      d=sqrt((IJ(:,1)-(i+di)).^2+(IJ(:,2)-(j-dj)).^2);
      iM=find(d==0);
      IJ(iM,:)=nan;
      in=find(~isnan(IJ(:,1)));
      IJ=IJ(in,:);

      chckp=0;
      if chckp>0
	clf;
	plot(IJ(:,1),IJ(:,2),'r.-');
	hold
	plot(i,j,'k*');
      end

  % Convert lon/lat -> cartesian km:
  % centered at i,j grid point:
      clear XY
      lon0=X(j,i);
      lat0=Y(j,i);
      for ik=1:length(IJ)
	i1=IJ(ik,1);
	j1=IJ(ik,2);
	lon1=X(j1,i1);
	lat1=Y(j1,i1);
  %      dd=distance_spheric_coord(lat1,lon1,lat0,lon0);
  %      alph=atan2((lat1-lat0),(lon1-lon0));
  %      dx=dd*cos(alph);
  %      dy=dd*sin(alph);
	dx=distance_spheric_coord(lat0,lon1,lat0,lon0);
%	sgnx=sign(lon1-lon0); <---- wrong for ASR flipped +lon direction in upper part of grid
        sgnx=sign(i1-i0);
	dy=distance_spheric_coord(lat1,lon0,lat0,lon0);
%	sgny=sign(lat1-lat0);
        sgny=sign(j1-j0);
	
	XY(ik,1)=sgnx*dx;
	XY(ik,2)=sgny*dy;
      end

      INDX(lI).IJ=IJ;
      INDX(lI).XY=XY;

    else
      IJ=INDX(lI).IJ;
      XY=INDX(lI).XY;
      if isempty(IJ), continue; end;
    end;
    nC=length(IJ)-1;
    
      
% Check how many nans on the contour:    
    clear uc vc
    for k=1:nC
      uc(k,1)=U(IJ(k,2),IJ(k,1));
      vc(k,1)=V(IJ(k,2),IJ(k,1));
    end;
    Nmax=round(0.2*nC);
    Inan=find(isnan(uc));
    if length(Inan)>Nmax, 
      CRC(j,i)=nan;
      continue; 
    end;

% To check:   
% Calculate vorticity in the domain in each point  
% should be somewhat close to zita= C/A (from circ. theorem), 
% at least signs
% should be similar
    chckv=logical(0);
    if chckv
      IN=inpolygon(IIND,JIND,IJ(:,1),IJ(:,2));
      II=find(IN==1);
      uu=U(II);
      vv=V(II);
      clear vrt
      
      for ipp=1:length(II);
	[ja,ia]=ind2sub([m,n],II(ipp));
	xa=X(ja,ia-1);
	ya=Y(ja,ia-1);
	xb=X(ja,ia+1);
	yb=Y(ja,ia+1);
	dx=distance_spheric_coord(ya,xa,yb,xb);
	dvdx=(V(ja,ia+1)-V(ja,ia-1))/dx;
	xa=X(ja-1,ia);
	ya=Y(ja-1,ia);
	xb=X(ja+1,ia);
	yb=Y(ja+1,ia);
	dy=distance_spheric_coord(ya,xa,yb,xb);
	dudy=(U(ja+1,ia)-U(ja-1,ia))/dy;
	vrt(ipp,1)=dvdx-dudy;
      end
      Cch=mean(vrt);
    end

    
%
% Calculate area of a planar non-self-intersecting polygon:
% http://mathworld.wolfram.com/PolygonArea.html
    A=0;
    for ik=1:nC    
      A=A+(XY(ik,1)*XY(ik+1,2)-XY(ik+1,1)*XY(ik,2));
    end
    A=0.5*abs(A);
    
% Circulation theorem:
    clear crc
    for ik=1:nC;
      i1=IJ(ik,1);
      j1=IJ(ik,2);
      i2=IJ(ik+1,1);
      j2=IJ(ik+1,2);
      x1=XY(ik,1);
      y1=XY(ik,2);
      x2=XY(ik+1,1);
      y2=XY(ik+1,2);
      
% Check orientation of the segment vector
% should go c/clockwise
      aa=[(x2-x1),(y2-y1),0];
      bb=[-x1,-y1,0]; % directed toward the center of the contour, x0=0,y0=0
      axb=cross(aa,bb);
      if axb(3)<0,
%	fprintf('Changing segment orientation j=%i, i=%i \n',j,i);
        i2=IJ(ik,1);
        j2=IJ(ik,2);
        i1=IJ(ik+1,1);
        j1=IJ(ik+1,2);
      end  
      
      u1=U(j1,i1);
      u2=U(j2,i2);
      v1=V(j1,i1);
      v2=V(j2,i2);
      x1=XY(ik,1);
      y1=XY(ik,2);
      x2=XY(ik+1,1);
      y2=XY(ik+1,2);
      crc(ik,1) = 0.5*((u2+u1)*(x2-x1)+(v2+v1)*(y2-y1));
    end

    C=nansum(crc);
    VRT(j,i)=C/A;  % s^-1
    AREA(j,i)=A;
%    if j==276 & i==206; keyboard; end;
%    keyboard
  end
end

VORT.VRT=VRT;
VORT.AREA=AREA;

return
