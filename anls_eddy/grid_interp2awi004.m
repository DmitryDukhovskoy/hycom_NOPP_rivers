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

YR1 = 2005;
YR2 = 2007;

rg=9806;  % convert pressure to depth, m
huge=1e10;
omg = 7.2921159e-5; 


sfig = 0;
zz0  = -100;  %  depth for interpolation

pthmat  = '/Net/ocean/ddmitry/HYCOM/ARCc/data_eddy_AWI/';
pthout  = '/Net/tholia/ddmitry/hycom/ARCc0.04/data_awi_intrp/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/topo_grid/';
pthfig  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/%s/fig_2D/',...
		  expt);

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


% Barycentric interpolation weights and HYCOM indices
fout    = sprintf('%sindx2awi004.mat',pthmat);
fprintf('Loading %s\n',fout);
load(fout);


ftopo = sprintf('%sdepth_ARCc0.04_%s.nc',pthtopo,TV); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

% Find AWI domain
dd=distance_spheric_coord(LAT,LON,lt1,ln1);
[j11,i11]=find(dd==min(min(dd)));
dd=distance_spheric_coord(LAT,LON,lt2,ln1);
[j21,i21]=find(dd==min(min(dd)));
dd=distance_spheric_coord(LAT,LON,lt1,ln2);
[j12,i12]=find(dd==min(min(dd)));
dd=distance_spheric_coord(LAT,LON,lt2,ln2);
[j22,i22]=find(dd==min(min(dd)));
JH1=min([j11,j12,j21,j22])-2;
JH2=max([j11,j12,j21,j22])+2;
IH1=min([i11,i12,i21,i22])-2;
IH2=max([i11,i12,i21,i22])+2;


f_map=0;
if f_map==1
  figure(1); clf;
  contour(HH,[0 0],'k'); hold
  contour(HH,[zz0 zz0],'r');
  contour(LON,[-170:10:170],'Color',[0.5 0.5 0.5]);
  contour(LAT,[40:10:79],'Color',[0.5 0.5 0.5]);
  axis('equal');
  lt1=min(alat);
  lt2=max(alat);
  ln1=min(alon);
  ln2=max(alon);

  contour(LON,[ln1 ln1],'g');
  contour(LON,[ln2 ln2],'g');
  contour(LAT,[lt1 lt1],'g');
  contour(LAT,[lt2 lt2],'g');
  title('AWI Fram domain');
  
      
end

f_chck=0;
if f_chck==1
  i0=100;
  j0=80;
  x0=alon(i0);
  y0=alat(j0);
  im1=INDX(1).hycom_I(j0,i0);
  jm1=INDX(1).hycom_J(j0,i0);
  im2=INDX(2).hycom_I(j0,i0);
  jm2=INDX(2).hycom_J(j0,i0);
  im3=INDX(3).hycom_I(j0,i0);
  jm3=INDX(3).hycom_J(j0,i0);
  
  lmb1=INDX(1).brctr_lmbd(j0,i0);
  lmb2=INDX(2).brctr_lmbd(j0,i0);
  lmb3=INDX(3).brctr_lmbd(j0,i0);
  
  figure(10); clf;
  plot(x0,y0,'r*');
  hold on;
  x1=LON(jm1,im1);
  x2=LON(jm2,im2);
  x3=LON(jm3,im3);
  y1=LAT(jm1,im1);
  y2=LAT(jm2,im2);
  y3=LAT(jm3,im3);
  plot([x1 x2],[y1 y2],'b-');
  plot([x1 x3],[y1 y3],'b-');
  plot([x3 x2],[y3 y2],'b-');
  text(x1,y1,sprintf('%4.3f',lmb1));
  text(x2,y2,sprintf('%4.3f',lmb2));
  text(x3,y3,sprintf('%4.3f',lmb3));
  
end

      
%    keyboard
% Interpolation into AWI grid:
im1 = INDX(1).hycom_I; % i HYCOM index of vertex 1
jm1 = INDX(1).hycom_J; % j HYCOM index of vertex 1
im2 = INDX(2).hycom_I;
jm2 = INDX(2).hycom_J;
im3 = INDX(3).hycom_I;
jm3 = INDX(3).hycom_J;
[ma,na] = size(im1);
im1 = reshape(im1,ma*na,1);
jm1 = reshape(jm1,ma*na,1);
im2 = reshape(im2,ma*na,1);
jm2 = reshape(jm2,ma*na,1);
im3 = reshape(im3,ma*na,1);
jm3 = reshape(jm3,ma*na,1);
I1  = sub2ind([mm,nn],jm1,im1);
I2  = sub2ind([mm,nn],jm1,im1);
I3  = sub2ind([mm,nn],jm1,im1);

lmb1 = INDX(1).brctr_lmbd;
lmb2 = INDX(2).brctr_lmbd;
lmb3 = INDX(3).brctr_lmbd;
lmb1 = reshape(lmb1,ma*na,1);
lmb2 = reshape(lmb2,ma*na,1);
lmb3 = reshape(lmb3,ma*na,1);

Hi = lmb1.*HH(I1)+lmb2.*HH(I2)+lmb3.*HH(I3);


Hi = reshape(Hi,[ma,na]);
Lmsk = Hi*0;
Lmsk(Hi<=zz0)=1;



f_chip = 0;
if f_chip==1
  pcolor(Hi); shading flat;
  colorbar
  axis('equal');
  contour(Hi,[0 0],'k')
  contour(Lmsk,[0.99 0.99],'b')
  contour(Hi,[-5000:500:-99],'Color',[0.4 0.4 0.4]);

  title('HYCOM Topo Interp  AWIgrid, blue - 100m');
%       caxis([-2 8]);
%  caxis([31 35]);
  set(gca,'xlim',[0 na],...
	  'ylim',[0 ma]);

  btx = 'grid_interp2awi004.m';
  bottom_text(btx,'pwd',1);

end


fintrp = sprintf('%sarc004_%s_grid2awi.mat',...
	       pthout,expt);
fprintf('Saving %s\n\n',fintrp);
save(fintrp,'Hi','Lmsk');


