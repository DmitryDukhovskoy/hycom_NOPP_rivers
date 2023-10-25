% Parallel code
% Plot CFSV2 - NCEP wind fields 0.3125 dgr res
%
% CFSv2 data are interpolated onto strethed
% bi-pole grid of HYCOM
% i.e., for a given point with coordinates p0,l0 on CFSv2 projection
% find corresponding i,j index on HYCOM grid
%
% Interpolation is carried out using barycentric coordinates for triangles
% Idea: For a pnt with coord. x0,y0 Find three closest points: x1,y1; 
% x2,y2; x3,y3; Get barycentric coordinates of the pnt x0,y0 and then get indices
% I,J for this point corresponding to HYCOM grid
% 
%


addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08/;
addpath /usr/people/ddmitry/codes/MyMatlab/;
startup;

clear
close all



% Parallel code:
delete(gcp('nocreate')); % make sure no parallel pool is running
pp = parpool(12);

year=2015;
mnth=12;
dd=1;
%pthdat=(['/Net/Movies0/ddmitry/vector_winds_ccmp_level30/',int2str(year),'/']);
pthd = '/Net/kronos/ddmitry/ncep_cfsv2/';
pth72='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.72/topo_grid/';
pth8='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';

%ffout=[pth8,'cfsv2_gridindx_arc008all.mat'];
%ffoutA=[pth8,'cfsv2_gridindx_arc008TMP.mat'];
ffout=[pth72,'cfsv2_gridindx_arc072all.mat'];
ffoutA=[pth72,'cfsv2_gridindx_arc072TMP.mat'];


date_s=sprintf('%4i%2.2i%2.2i',year,mnth,dd)
fnm = sprintf('cfsv2-sec2_%4.4i_mon_uv-10m.nc',year);
fp = [pthd,fnm];

% Flags:
f_topo=1;       % = 1 - use ARCc0.72 model topography,
                % = 2 - use ETOPO-2 topography
                % = 3 - use ARCc0.08 topography
f_smat = 1;	% =1 - save final mat file; =2 save temporary and
                % final; =0 - do not save mat files   
f_load_mat=0;   %=1- load previous indx, start from last record 

if f_smat==0
  fprintf('output mat files will not be created\n\n');
end


elon = nc_varget(fp,'Longitude');
alat = nc_varget(fp,'Latitude');
n=length(elon);
m=length(alat);

% Region of interest:
% Longitudes in the data file are from 0.125 through 359.875
% Nordic seas:
%y1=55;
%y2=85;
%x1=360-45.125;
%x2=70;
% Whole AO domain
% tricky for the ARCc0.08 grid:
% easier to specify last/first
% index in ARCc grid:
y1=9999;
y2=9999;
x1=9999;
x2=9999;


dx=abs(elon(2)-elon(1));

if y1<1e3
  dst=sqrt((elon-x1).^2);
  i1=min(find(dst==min(dst)));
  dst=sqrt((alat-y1).^2);
  j1=min(find(dst==min(dst)));
  dst=sqrt((elon-x2).^2);
  i2=min(find(dst==min(dst)));
  dst=sqrt((alat-y2).^2);
  j2=min(find(dst==min(dst)));
else
  i1=1;
  i2=length(elon);
  dst=sqrt((alat-60).^2);
  j1=min(find(dst==min(dst)));
  j2=length(alat);
end


if x1>x2   % the region is cut by 0 meridian
  e1=elon(i1:end);
  e2=elon(1:i2);
  X=[e1;e2];
else
  X=elon(i1:i2);
end;
Y=alat(j1:j2);

X(X>180)=X(X>180)-360;


% Topography:

% Get etopo2 data:

  xs=min(X);
  xe=max(X);
  ys=y1;
  ye=y2;

if f_topo==2
  disp('Extracting etopo2 bath for the Arctic ...');
  [eltopo,altopo,mtopo,ntopo,H]=get_grid;
  dp=4;
elseif f_topo==1
  disp('Extracting topography for ARCc0.72');
  fsv=[pth72,'new_bath072.mat'];
  load(fsv);
  eltopo=elon;
  altopo=alat;
  H=hnew;
  H(isnan(H))=10;
  clear hnew elon alat;
  dp=1;
elseif f_topo==3
  disp('Extracting topography for ARCc0.08');
  fsv=[pth8,'bath_arc008.mat'];
  load(fsv);
  eltopo=plon;
  altopo=plat;
  H=h;
  clear h plon plat
  dp=1;
end

[at1,at2]=size(H);
[mm,nn] = size(H);
% Extract bathymetry for the region:
% Note that ETOPO projection is Polar while
% CCMP data are Lambert
%  dst=sqrt((eltopo-xs).^2+(altopo-ys).^2);
if x1<1e3
dst=distance_spheric_coord(ys,xs,altopo,eltopo);
[jtp1,itp1]=find(dst==min(min(dst)));
%  dst=sqrt((eltopo-xe).^2+(altopo-ys).^2);
dst=distance_spheric_coord(ye,xe,altopo,eltopo);
[jtp2,itp2]=find(dst==min(min(dst)));
else % take all topo domain
  jtp1=1;
  jtp2=mm;
  itp1=1;
  itp2=nn;
end


XT=eltopo(jtp1:dp:jtp2,itp1:dp:itp2);
YT=altopo(jtp1:dp:jtp2,itp1:dp:itp2);
HT=H(jtp1:dp:jtp2,itp1:dp:itp2);
clear eltopo altopo H

[mtopo,ntopo]=size(XT);
[IT,JT]=meshgrid((1:ntopo),(1:mtopo));

%[XX,YY]=meshgrid(X,Y);

ii1=1;
if f_load_mat>0
  fprintf('Loading %s\n',ffoutA);
  fprintf('Continue indices search\n');
  load(ffoutA);
  ii1=ii+1;
  fprintf('Starting from ii=%i\n',ii1);
end


% Find indeces of TOPO grid corresponding to CFSv2 wind data:
% Project CFSv2 projection onto Polar grid of Topo data
% ii=574, jj=21 or 66
disp('Matching indices ...');
nw=length(X);
mw=length(Y);
for ii=ii1:nw
%for ii=131:131
%  disp(['Searching indices, ',int2str(ii)]);
  if mod(ii,4)==0;
    rat=ii/nw*100;
    fprintf('Searching indices %5.1f %%...\n\n',rat);
  end

  
% Parallel loop:  
  clear iwp jwp
%  for jj=4:4
  parfor jj=1:mw
%    fprintf('jj=%i\n',jj);
%    disp(jj);
    x0=X(ii);
    y0=Y(jj);
% Check if pnt is inside hycom/topo domain
% CFSv2 points on the boundart of
% the hycom/topo grid are discarded
% to avoid possible errors
% in finding 3 surrounding points
% First check if the point is completely outside:
    Ix=find(XT>=x0-0.15 & XT<=x0+0.15);
     if isempty(Ix)
      iwp(jj,1)=-x0;
      jwp(jj,1)=-y0;
      continue;
    end
    Ymn=min(YT(Ix));
    if y0<min(Ymn)
      iwp(jj,1)=-x0;
      jwp(jj,1)=-y0;
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
      iwp(jj,1)=-x0;
      jwp(jj,1)=-y0;
      continue;
    end      
    


    
%    fprintf('%i, %i x0=%5.2f y0=%5.2f\n',jj,ii,x0,y0);
    xm1=[];
    ym1=[];
    xm2=[];
    ym2=[];
    xm3=[];
    ym3=[];
%    clear XYT
%    dst=sqrt((XT-x0).^2+(YT-y0).^2);
    dst=distance_spheric_coord(YT,XT,y0,x0);
    [j0,i0]=find(dst==min(min(dst)));


    j1=j0;
    i1=i0;
    xm1=XT(j1,i1);
    ym1=YT(j1,i1);

%      dst=sqrt((XT-x0).^2+(YT-y0).^2);
    dst=distance_spheric_coord(YT,XT,y0,x0);
    dst(j1,i1)=1e6;
    [j2,i2]=find(dst==min(min(dst)));
    dst(j2,i2)=1e6;

%   3rd vertex can not be on the same line with vertices 1, 2
    if (j2==j1)
      dst(j1,:)=1e6;
    elseif (i2==i1)
      dst(:,i1)=1e6;
    end

    [j3,i3]=find(dst==min(min(dst)));
    xm2=XT(j2,i2);
    ym2=YT(j2,i2);
    xm3=XT(j3,i3);
    ym3=YT(j3,i3);

%      XYT(1,1)=xm1;
%      XYT(1,2)=ym1;
%      XYT(2,1)=xm2;
%      XYT(2,2)=ym2;
%      XYT(3,1)=xm3;
%      XYT(3,2)=ym3;
    XYT=[xm1,ym1;xm2,ym2;xm3,ym3];

%      fprintf('1: jj=%i ii=%i \n',jj,ii);
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
      [j3,i3]=find(dst==min(min(dst)));
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
      

    if Ih>ntopo | Jh>mtopo
      disp(['Pnt from ii=',int2str(ii),' jj=',int2str(jj)]);
      disp(['Assigned i=',int2str(Ih),' should be <  ',int2str(ntopo)]);
      disp(['Assigned j=',int2str(Jh),' should be <  ',int2str(ntopo)]);
      error('*** ERR: Assigned index is outside HYCOM grid');
   end

%      IWP(jj,ii)=Ih;
%      JWP(jj,ii)=Jh;
     iwp(jj,1)=Ih;
     jwp(jj,1)=Jh;

    plt_tmp=1;
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

  niwp=length(iwp);
  IWP(1:niwp,ii)=iwp;
  JWP(1:niwp,ii)=jwp;

  if f_smat>1
    fprintf('Saving tmpr: %s\n',ffoutA);
    save(ffoutA,'ii','*WP','HT','XT','YT','X','Y');   
  end
  
  
end;    % for ii

%matlabpool close;

delete(gcp('nocreate')); % make sure no parallel pool is running


disp('All indices have been identified ...');
%pause
  
%
% IWP,JWP - indices for wind vectors on topo grid
% HT - topography
% XT, YT - lon/lat for topography
% X, Y - coordinates of selected region of wind data
%disp(['Saving interpolated grid and topography ...']);
%ffout=[pth8,'cfsv2_gridindx_arc008.mat'];
%ffout=[pth8,'cfsv2_gridindx_arc008all.mat'];
if f_smat>0
  fprintf('Saving %s\n',ffout);
  save(ffout,'ii','*WP','HT','XT','YT','X','Y');   
end


% ================================





