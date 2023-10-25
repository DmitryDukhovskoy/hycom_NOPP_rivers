function VRT = sub_vorticity(X,Y,U,V,RR,INDX);
%
% Calculate vorticity
% by standard approach 
% using definition of the vorticity:
%  zt = dv/dx - du/dy
% then vorticity is averaged
% over an area of a rectangular region 
% with half-side= RR
%
% Commented out: Also calculate ~circulation around the region
% to see if it closed or not
% with 0 - indicating circular flow around a point (closed)
%          i.e. flows across the point across X and Y
%          compensate each other (scaled to be +1 or -1)
% 
%         <-  -1
%                  ^
%      |    *(i,j) | +1
%      V
%    -1     ->  +1
%
% Input:
%   U,V - wind field
%   X,Y - geogr. coordinates of U,V
%         Note that U,V, X, Y has to be
%         "glued" along the 0 meridian
%         so that X(1)<X(2)... <X(N), where
%         X(1) is the westernmost location
%   RR - averaged ared size ("diameter"), km, will be rounded 
%        to the nearest grid points
%
% Output: VORT - structured array
%

% Convert RR to m if needed:
if RR<1e3
  RR=RR*1e3;
end

l=min(size(X));
if l==1
  [X,Y]=meshgrid(X,Y);
end;

[m,n]=size(U);

DX=zeros(m,n);
DY=zeros(m,n);
ln1=X(:,1:n-1);
ln2=X(:,2:n);
lt1=Y(:,1:n-1);
lt2=Y(:,2:n);
DX0=distance_spheric_coord(lt1,ln1,lt2,ln2); % m
DX0(:,n)=DX0(:,n-1);

ln1=X(1:m-1,:);
ln2=X(2:m,:);
lt1=Y(1:m-1,:);
lt2=Y(2:m,:);
DY0=distance_spheric_coord(lt1,ln1,lt2,ln2); % m
DY0(m,:)=DY0(m-1,:);

% Find # of grid points for averaging
DD=sqrt(DX0.^2+DY0.^2);
Smn=nanmean(nanmean(DD));
NP=round(RR/Smn);
if NP==0; NP=1; end; % RR< grid size

% 

VRT=zeros(m,n)*nan;
DU=VRT;
DV=VRT;
%keyboard
% Distance between differenced points:
DX=zeros(m,n)*nan;
DY=zeros(m,n)*nan;
DU=DX;
DV=DX;
ln1=X(:,1:n-2*NP);
ln2=X(:,1+2*NP:n);
lt1=Y(:,1:n-2*NP);
lt2=Y(:,1+2*NP:n);
dmm=distance_spheric_coord(lt1,ln1,lt2,ln2); % m
DX(:,1+NP:n-NP)=dmm;

ln1=X(1:m-2*NP,:);
ln2=X(1+2*NP:m,:);
lt1=Y(1:m-2*NP,:);
lt2=Y(1+2*NP:m,:);
dmm=distance_spheric_coord(lt1,ln1,lt2,ln2); % m
DY(1+NP:m-NP,:)=dmm;

uj1=U(1:m-2*NP,:);
uj2=U(1+2*NP:m,:);
vi1=V(:,1:n-2*NP);
vi2=V(:,1+2*NP:n);
dmm=uj2-uj1;
DU(1+NP:m-NP,:)=dmm;
dmm=vi2-vi1;
DV(:,1+NP:n-NP)=dmm;


DUDY=DU./DY;
DVDX=DV./DX;
VRT = DVDX-DUDY;
crY=VRT*nan;
crX=crY;
% Calculate "circulation" 
%dmm=uj1+uj2;
%crY(1+NP:m-NP,:)=dmm;
%dmm=vi1+vi2;
%crX(:,1+NP:n-NP)=dmm;
%CRC=crX+crY;

%keyboard

return