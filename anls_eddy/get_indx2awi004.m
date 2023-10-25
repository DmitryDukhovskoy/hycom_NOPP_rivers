% Find interpolation indices 
% for interpolating ARCc0.05 -> AWI
% Using barycentirc interpolation:
% p(x,y)=lmb1*f(x1,y1)+lmb2*f(x2,y2)+lmb3*f(x3,y3)
% Interpolate HYCOM 0.04 fields onto
% AWI grid
% for calculating vorticity etc
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

expt = '011';
TV = '17DD';
rg=9806;  % convert pressure to depth, m
huge=1e10;
omg = 7.2921159e-5; 
dnmb = datenum(2005,8,21); % date to plot

sfig = 0;
zz0  = -100;  %  depth for interpolation

pthmat  = '/Net/ocean/ddmitry/HYCOM/ARCc/data_eddy_AWI/';
pthout  = '/Net/tholia/ddmitry/hycom/ARCc0.04/data_awi_intrp/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/topo_grid/';
pthfig  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/%s/fig_2D/',...
		  expt);
fout    = sprintf('%sindx2awi004.mat',pthmat);


% Get AWI grid:
fmat = sprintf('%sgrid_eddy_tracking.mat',pthmat);
fprintf('Loading %s\n',fmat);
AWI = load(fmat);
alat=AWI.lat;
alon=AWI.lon;
lt1=min(alat);
lt2=max(alat);
ln1=min(alon);
ln2=max(alon);


ftopo = sprintf('%sdepth_ARCc0.04_%s.nc',pthtopo,TV); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

f_map=0;
if f_map==1
figure(1); clf;
  contour(HH,[0 0],'k'); hold
  contour(HH,[zz0 zz0],'r');
  contour(LON,[-170:10:170],'Color',[0.5 0.5 0.5]);
  contour(LAT,[40:10:79],'Color',[0.5 0.5 0.5]);
  axis('equal');

  contour(LON,[ln1 ln1],'g');
  contour(LON,[ln2 ln2],'g');
  contour(LAT,[lt1 lt1],'g');
  contour(LAT,[lt2 lt2],'g');
  title('AWI Fram domain');
end

% subsample smaller domain
dd=distance_spheric_coord(LAT,LON,lt1,ln1);
[j11,i11]=find(dd==min(min(dd)));
dd=distance_spheric_coord(LAT,LON,lt2,ln1);
[j21,i21]=find(dd==min(min(dd)));
dd=distance_spheric_coord(LAT,LON,lt1,ln2);
[j12,i12]=find(dd==min(min(dd)));
dd=distance_spheric_coord(LAT,LON,lt2,ln2);
[j22,i22]=find(dd==min(min(dd)));
J1=min([j11,j12,j21,j22])-3;
J2=max([j11,j12,j21,j22])+3;
I1=min([i11,i12,i21,i22])-3;
I2=max([i11,i12,i21,i22])+3;


YT=LAT(J1:J2,I1:I2);
XT=LON(J1:J2,I1:I2);
HS=HH(J1:J2,I1:I2);

fprintf('Start seraching indices ...\n');

na=length(alon);
ma=length(alat);
nma=na*ma;
cnc=0;
tic;
for ii=1:na
  for jj=1:ma
    cnc=cnc+1;
    if mod(cnc,1000)==0
      fprintf('Finding indices, %6.2f%% done %6.2f min...\n',...
	      cnc/nma*100,toc/60);
      tic;
    end
    
    x0=alon(ii);
    y0=alat(jj);
% Find 4 closest points in the verticies
    xm1=[];
    ym1=[];
    xm2=[];
    ym2=[];
    xm3=[];
    ym3=[];

    dst=distance_spheric_coord(YT,XT,y0,x0);
    [j1,i1]=find(dst==min(min(dst)));
    xm1=XT(j1,i1);
    ym1=YT(j1,i1);
    dst(j1,i1)=1e9;
    
    [j2,i2]=find(dst==min(min(dst)));
    xm2=XT(j2,i2);
    ym2=YT(j2,i2);
    dst(j2,i2)=1e9;
    
    [j3,i3]=find(dst==min(min(dst)));
    xm3=XT(j3,i3);
    ym3=YT(j3,i3);
    dst(j3,i3)=1e9;
%
% Calculate distances instead of geogr. coord
%    xd1=0;
%    yd1=0;
%    xd2=distance_spheric_coord(ym1,xm1,ym1,xm2);
%    yd2=distance_spheric_coord(ym1,xm1,ym2,xm1);
%    xd3=distance_spheric_coord(ym1,xm1,ym1,xm3);
%    yd3=distance_spheric_coord(ym1,xm1,ym3,xm1);
%    xd0=distance_spheric_coord(ym1,xm1,ym1,x0);
%    yd0=distance_spheric_coord(ym1,xm1,y0,xm1);
    
    XYT=[xm1,ym1;xm2,ym2;xm3,ym3];
    [lmb1,lmb2,lmb3] = barycentric_coord(XYT,x0,y0);
%    XYT2=[xd1,yd1;xd2,yd2;xd3,yd3];
%    [lb1,lb2,lb3] = barycentric_coord(XYT2,xd0,yd0); % very close to above one
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
%	xm3=-999;
%	ym3=-999;
%        break;
	error('*** ERR: stopping ...');
      end
      dst(j3,i3)=1e6;
      [j3,i3]=find(dst==min(min(dst)));
      xm3=XT(j3,i3);
      ym3=YT(j3,i3);
      XYT(3,1)=xm3;
      XYT(3,2)=ym3;
%      plot(xm3,ym3,'m+');
      [lmb1,lmb2,lmb3] = barycentric_coord(XYT,x0,y0);
      if lmb1>=0 & lmb1<=1 & lmb2>=0 & lmb2<=1 & lmb3>=0 & lmb3<=1
	LMB=0;
      end;
    end
    
    plt_tmp=0;
    if plt_tmp==1
      figure(10); clf;
      plot(XT,YT,'c.');
      axis('equal');
      hold on
      plot(x0,y0,'r*');
      plot(xm1,ym1,'ro');
      plot(xm2,ym2,'go');
      plot(xm3,ym3,'bo');
    end;
    
    INDX(1).hycom_I(jj,ii)    = i1+I1-1;
    INDX(1).hycom_J(jj,ii)    = j1+J1-1;
    INDX(1).brctr_lmbd(jj,ii) = lmb1;

    INDX(2).hycom_I(jj,ii)    = i2+I1-1;
    INDX(2).hycom_J(jj,ii)    = j2+J1-1;
    INDX(2).brctr_lmbd(jj,ii) = lmb2;

    INDX(3).hycom_I(jj,ii)    = i3+I1-1;
    INDX(3).hycom_J(jj,ii)    = j3+J1-1;
    INDX(3).brctr_lmbd(jj,ii) = lmb3;
    
  end
end

fprintf('Saving %s\n',fout);
save(fout,'INDX');

