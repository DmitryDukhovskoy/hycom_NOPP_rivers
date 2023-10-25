% EN4 data are interpolated onto 
% bi-polar grid of HYCOM
% i.e., for a given point with coordinates p0,l0 on EN4 projection
% find corresponding i,j index on HYCOM grid
%
% Interpolation is carried out using barycentric coordinates for triangles
% Idea: For a pnt with coord. x0,y0 Find three closest HYCOM points: x1,y1; 
% x2,y2; x3,y3; 
% Get barycentric coordinates of the pnt x0,y0 and then get HYCOM indices
% I,J for this EN4 point 
% 
%


addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08/;
addpath /usr/people/ddmitry/codes/MyMatlab/;
startup;

clear
close all

f_smat = 1;

pthmat  = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/data_mat/';
pthin   = '/Net/data/GDEM4/';  % climatology data with no land
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
ffout   = sprintf('%sGDEM_gridindx_arc008.mat',pthmat);

ftopo = sprintf('%sdepth_ARCc0.08_07.nc',pthtopo); 
HT  = nc_varget(ftopo,'Bathymetry'); 
XT = nc_varget(ftopo,'Longitude');
YT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(HT); 
[IT,JT] = meshgrid((1:nn),(1:mm));

[DX,DY]=sub_dx_dy(XT,YT);


%-------------------------------
im=3;
finT=sprintf('%sptgdemv4f%2.2i.nc4',pthin,im);  % ptgdemv4f##.nc4
ZZ = nc_varget(finT,'Depth');
TT = nc_varget(finT,'Potential_Temperature');
LON = nc_varget(finT,'Longitude');
LAT = nc_varget(finT,'Latitude');

I = find(LON>180);
LON(I)=LON(I)-360;

[X,Y] = meshgrid(LON,LAT);
[mw,nw] = size(X);
Ib = find(Y>37);
nI = length(Ib);

for ik=1:nI
  i0 = Ib(ik);
  [jj,ii]=ind2sub(size(X),i0);

  x0=X(jj,ii);
  y0=Y(jj,ii);
  
% Check if pnt is inside hycom/topo domain
% CFSR points on the boundart of
% the hycom/topo grid are discarded
% to avoid possible errors
% in finding 3 surrounding points
% First check if the point is completely outside:
  Ix=find(XT>=x0-0.15 & XT<=x0+0.15);
  if isempty(Ix)
    IWP(jj,ii)=-x0;
    JWP(jj,ii)=-y0;
    continue;
  end
  Ymn=min(YT(Ix));
  if y0<min(Ymn)
    IWP(jj,ii)=-x0;
    JWP(jj,ii)=-y0;
    continue;
  end
  a=XT(Ix);
  b=YT(Ix);

% Second, Corner/bndry  points - discard    
  I1=find(a<=x0);
  I2=find(a==x0);
  I3=find(a>x0);
 b=YT(Ix);
  J1=find(b<=y0);
  J2=find(b==y0);
  J3=find(b>y0);

  if isempty(I1) | isempty(I3) | isempty(J1) | isempty(J3)
    IWP(jj,ii)=-x0;
    JWP(jj,ii)=-y0;
    continue;
  end      

  xm1=[];
  ym1=[];
  xm2=[];
  ym2=[];
  xm3=[];
  ym3=[];
%    clear XYT
%    dst=sqrt((XT-x0).^2+(YT-y0).^2);
  dst=distance_spheric_coord(YT,XT,y0,x0);
  [j0,i0]=find(dst==min(min(dst)),1);

  j1=j0;
  i1=i0;
  xm1=XT(j1,i1);
  ym1=YT(j1,i1);

%      dst=sqrt((XT-x0).^2+(YT-y0).^2);
  dst=distance_spheric_coord(YT,XT,y0,x0);
  dst(j1,i1)=1e6;
  [j2,i2]=find(dst==min(min(dst)),1);
  dst(j2,i2)=1e6;

%   3rd vertex can not be on the same line with vertices 1, 2
  if (j2==j1)
    dst(j1,:)=1e6;
  elseif (i2==i1)
    dst(:,i1)=1e6;
  end

  [j3,i3]=find(dst==min(min(dst)),1);
  xm2=XT(j2,i2);
  ym2=YT(j2,i2);
  xm3=XT(j3,i3);
  ym3=YT(j3,i3);

  XYT=[xm1,ym1;xm2,ym2;xm3,ym3];

  [lmb1,lmb2,lmb3] = barycentric_coord(XYT,x0,y0);
%
% Check if the pnt is inside the triangle
% if not take another points
% all lambdas has to be >0 and <1 
  LMB=1;
  if lmb1>=0 & lmb1<=1 & lmb2>=0 & lmb2<=1 & lmb3>=0 & lmb3<=1
    LMB=0;
  end;
  cntr=0;
  while abs(LMB)>0;
    cntr=cntr+1;
    if cntr>20
      fprintf('WARNING: ii=%i, jj=%i\n',ii,jj);
      fprintf('WARNING: could not locate 3rd vortex\n');
      fprintf('%i, %i x0=%5.2f y0=%5.2f\n',jj,ii,x0,y0);
      fprintf('The point is outside the domain');
      fprintf('LMB=%d\n\n',LMB);
      xm3=-999;
      ym3=-999;
      break;
%	error('*** ERR: stopping ...');
    end
    dst(j3,i3)=1e6;
    [j3,i3]=find(dst==min(min(dst)),1);
    xm3=XT(j3,i3);
    ym3=YT(j3,i3);
    XYT(3,1)=xm3;
    XYT(3,2)=ym3;
%      plot(xm3,ym3,'m+');

%        fprintf('2: jj=%i ii=%i \n',jj,ii);
    [lmb1,lmb2,lmb3] = barycentric_coord(XYT,x0,y0);
    if lmb1>=0 & lmb1<=1 & lmb2>=0 & lmb2<=1 & lmb3>=0 & lmb3<=1
      LMB=0;
    end;
%        fprintf('=====      \n',jj,ii);

  end
  fprintf('%i, %i x0=%5.2f y0=%5.2f\n',jj,ii,x0,y0);
  fprintf('%i, %i xm3=%5.2f ym3=%5.2f\n',jj,ii,xm3,ym3);
% Indices of x0,y0 in HYCOM grid:
  if abs(LMB>0), % pnt is outside domain
    Ih=-x0;
    Jh=-y0;
  else
    Ih=lmb1*i1+lmb2*i2+lmb3*i3;
    Jh=lmb1*j1+lmb2*j2+lmb3*j3;
  end


  if Ih>nn | Jh>mm
    disp(['Pnt from ii=',int2str(ii),' jj=',int2str(jj)]);
    disp(['Assigned i=',int2str(Ih),' should be <  ',int2str(ntopo)]);
    disp(['Assigned j=',int2str(Jh),' should be <  ',int2str(ntopo)]);
    error('*** ERR: Assigned index is outside HYCOM grid');
 end

    IWP(jj,ii)=Ih;
    JWP(jj,ii)=Jh;
%     iwp(jj,1)=Ih;
%     jwp(jj,1)=Jh;

  plt_tmp=0;
  if plt_tmp==1
    plot(XT,YT,'b.');
    axis('equal');
    hold on
    plot(x0,y0,'r*');
    plot(xm1,ym1,'ro');
    plot(xm2,ym2,'go');
    plot(xm3,ym3,'mo');
  end;

%

%    keyboard

end;   % for jj

  
INDX.source = 'hycom_NOPP_rivers/EN4_anls/getHYCOMindices_GDEM.m';
INDX.IWP  = IWP;
INDX.JWP  = JWP;
INDX.XT   = XT;
INDX.YT   = YT;
INDX.HT   = HT;

if f_smat>0
  fprintf('Saving tmpr: %s\n',ffout);
  save(ffout,'INDX');   
end
