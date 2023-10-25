% GLBb0.08 19.0/19.1 reanalysis uses T07
% Arctic regional HYCOM/CICE uses T09:
% ------------------------------------------
% The difference between GLBb0.08 topography version 07 and 08 
% is that the latter was made to be compatible with CICE. 
% Thus, as you say, there is no regional.cice for version 07.
% Then the difference between 08 and 09 were 
% some changes to the Yenisei River in northern Russia. 
% It looks like we changed some ocean points in 08 to land points in 09.
% Joe
% ------------------------------------------
% Compare topo difference in the relax zone
% of ARCc 
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close

PTH.topo='/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
fltopo=sprintf('%sdepth_ARCc0.08_09.nc',PTH.topo);
HH   = nc_varget(fltopo,'Bathymetry');
alat = nc_varget(fltopo,'Latitude');
elon = nc_varget(fltopo,'Longitude');
LAT  = alat;
LON  = elon;
[m,n]= size(HH);
[mm,nn]= size(HH);

% Get T07:
% Read regional grid:
pthglb = '/nexsan/GLBb0.08/GLBb0.08_190/topo/';
frga   = sprintf('%sregional.grid.a',pthglb);
frgb   = sprintf('%sregional.grid.b',pthglb);
fdptha = sprintf('%sdepth_GLBb0.08_07.a',pthglb);
fdpthb = sprintf('%sdepth_GLBb0.08_07.b',pthglb);

fidRGb = fopen(frgb,'r');  % read I,J from regional.grid.b
aa  = fgetl(fidRGb);
dmm = aa(2:8);
IDM = str2num(dmm);
aa = fgetl(fidRGb);
dmm = aa(2:8);
JDM = str2num(dmm);
IJDM = IDM*JDM;
fclose(fidRGb);
npad=4096-mod(IJDM,4096);

fprintf('IDM=%i, JDM=%i\n',IDM,JDM);

% read lon/lat from GLBb regional grid file
fidRGa = fopen(frga,'r');
[plon,count] = fread(fidRGa,IJDM,'float32','ieee-be');
fseek(fidRGa,4*(npad+IJDM),-1);
[plat,count] = fread(fidRGa,IJDM,'float32','ieee-be');

disp('Reading lat/lon for GLBb0.08 ...')
plon=(reshape(plon,IDM,JDM))';
plat=(reshape(plat,IDM,JDM))';

fclose(fidRGa);

[mg,ng]=size(plon);

I=find(plon>180);
plon(I)=plon(I)-360;



% Find Arctic grid in GLBb grid:
D=distance_spheric_coord(LAT(1,1),LON(1,1),plat,plon);
[j1,i1]=find(D==min(min(D)));
fprintf('Closest point in GLBb is %8.4f m away from ARCc(1,1)\n',min(min(D)));
D=distance_spheric_coord(LAT(1,end),LON(1,end),plat,plon);
[j2,i2]=find(D==min(min(D)));
fprintf('Closest point in GLBb is %8.4f m away from ARCc(1,end)\n',min(min(D)));
D=distance_spheric_coord(LAT(end,end),LON(end,end),plat,plon);
[j3,i3]=find(D==min(min(D)));
fprintf('Closest point in GLBb is %8.4f m away from ARCc(end,end)\n',min(min(D)));
D=distance_spheric_coord(LAT(end,1),LON(end,1),plat,plon);
[j4,i4]=find(D==min(min(D)));
fprintf('Closest point in GLBb is %8.4f m away from ARCc(end,1)\n',min(min(D)));


% read bathymetry from regional.depth.a
depth_fid=fopen(fdptha,'r');

%dmm=fread(depth_fid,6,'float32','ieee-be');
%aa=fread(depth_fid,IJDM,'float32','ieee-be');
%dmm=fread(depth_fid,npad,'float32','ieee-be');
%fseek(depth_fid,6*4*(npad+IJDM),-1) % <--- not needed
[h,count]=fread(depth_fid,IJDM,'float32','ieee-be');
y=find(h>1e10);
h(y)=nan;
h=reshape(h,IDM,JDM)';

h=-h;
h(isnan(h))=100;
HG=h; % GLBb topo

fclose(depth_fid);

% Cut and rotate GLBb grid to match ARCc
figure(1);
contour(h,[0 0],'k');
hold;
plot(i1,j1,'r*');
plot(i2,j2,'b*');
plot(i3,j3,'m*');
plot(i4,j4,'g*');

pltb=0;
if pltb>0
  jj=1;
  for ii=100:100:n
    fprintf('i=%i\n',ii);
    D=distance_spheric_coord(LAT(jj,ii),LON(jj,ii),plat,plon);
    [j0,i0]=find(D==min(min(D)));
    plot(i0,j0,'g.');
  end;

  jj=mm;
  for ii=100:100:n
    fprintf('i=%i\n',ii);
    D=distance_spheric_coord(LAT(jj,ii),LON(jj,ii),plat,plon);
    [j0,i0]=find(D==min(min(D)));
    plot(i0,j0,'g.');
  end;

  ii=1;
  for jj=100:100:m
    D=distance_spheric_coord(LAT(jj,ii),LON(jj,ii),plat,plon);
    [j0,i0]=find(D==min(min(D)));
    fprintf('ii=%i, GLBb j0=%i, i0=%i\n',ii,j0,i0);
    plot(i0,j0,'g.');
  end;
  
  ii=nn;
  for jj=100:100:m
    D=distance_spheric_coord(LAT(jj,ii),LON(jj,ii),plat,plon);
    [j0,i0]=find(D==min(min(D)));
    fprintf('ii=%i, GLBb j0=%i, i0=%i\n',ii,j0,i0);
    plot(i0,j0,'g.');
  end;
end

axis('equal');
set(gca,'ylim',[1500 mg]);
dm1=plon;
dm1(dm1>175)=nan;
dm1(dm1<-175)=nan;
dm1(plat>89.)=nan;
contour(dm1,[-170:10:170],'y');
contour(dm1,[0 0],'r');
contour(plat,[35:5:89],'c');
contour(plat,[60 60],'r');
ipv=i3+nn/2-1;  % middle point in ARCc grid:
plot([ipv ipv],[j3 mg],'m--');



% the last row in GLBb - where ARCc is glued:
icup=ceil(0.5*(i1+i2)); % central pnt
icbt=ceil(0.5*(i3+i4));
% 
% Section that does not rotate/swap:
latup=plat(mg,i1:i2); % upper row in section that does not rotate
lonup=plon(mg,i1:i2);
latup1=plat(mg-1,i1:i2); % upper row 
lonup1=plon(mg-1,i1:i2);
% Section that is rotated and swapped:
% its top rows become bottom after rotation
% 2 rows overalp with the top 2 rows in the other section:
latbt=plat(mg,i3:i4); % btm row, topmost
lonbt=plon(mg,i3:i4);
latbt1=plat(mg-1,i3:i4); % btm row, next to top
lonbt1=plon(mg-1,i3:i4);
%lonbt2=fliplr(lonbt2);
% i1,j1 - lower left corner on ARCc grid
%

figure(3); clf;
contour(HH,[0 0],'k');
hold on;
dm1=LON;
dm1(dm1>175)=nan;
dm1(dm1<-175)=nan;
dm1(LAT>89.)=nan;
contour(dm1,[-170:10:170],'y');
contour(dm1,[0 0],'r');
contour(LAT,[35:5:89],'c');
contour(LAT,[60 60],'r');

ppnt=0;
if ppnt>0
  for ii=1:20:nn
    D=distance_spheric_coord(latup(ii),lonup(ii),LAT,LON);
    [j0,i0]=find(D==min(min(D)));
    fprintf('ii=%i, ARCc j0=%i, i0=%i\n',ii,j0,i0);
    plot(i0,j0,'b.');
  end;

  for ii=1:20:nn
    D=distance_spheric_coord(latbt(ii),lonbt(ii),LAT,LON);
    [j0,i0]=find(D==min(min(D)));
    fprintf('ii=%i, ARCc j0=%i, i0=%i\n',ii,j0,i0);
    plot(i0,j0,'g.');
  end;
end

% In GLBb, when pivoting 2 top rows are copied:
% segment1: goes to bottom after pivoting the GLBb grid
%           goes south to north along longitude ~79E
%           indices: i3:icup, 
% segment2: goes to bottom after pivoting
%           north-south along long ~ -105W
% segment3: Top after pivoting (does not change)
%           south-north along ~-105E
% segment4: Top does not pivot
%           north-south along 79E
% Match top and bottom segments on GLBb grid
sgm1_lon_r2 = fliplr(lonbt(1:nn/2));
sgm1_lat_r2 = fliplr(latbt(1:nn/2));
sgm1_lon_r1 = fliplr(lonbt1(1:nn/2));
sgm1_lat_r1 = fliplr(latbt1(1:nn/2));

sgm2_lon_r2 = fliplr(lonbt(nn/2+1:end));
sgm2_lat_r2 = fliplr(latbt(nn/2+1:end));
sgm2_lon_r1 = fliplr(lonbt1(nn/2+1:end));
sgm2_lat_r1 = fliplr(latbt1(nn/2+1:end));

sgm3_lon_r2 = lonup(1:nn/2);
sgm3_lat_r2 = latup(1:nn/2);
sgm3_lon_r1 = lonup1(1:nn/2);
sgm3_lat_r1 = latup1(1:nn/2);

sgm4_lon_r2 = lonup(nn/2+1:end);
sgm4_lat_r2 = latup(nn/2+1:end);
sgm4_lon_r1 = lonup1(nn/2+1:end);
sgm4_lat_r1 = latup1(nn/2+1:end);

figure(4); clf;
subplot(2,2,1);
hold on;
plot(sgm1_lon_r2,sgm1_lat_r2,'b.');
plot(sgm4_lon_r1,sgm4_lat_r1,'ro');
title('Segment 1-row2  and Segment 4-row 1'); 

subplot(2,2,2);
hold on;
plot(sgm1_lon_r1,sgm1_lat_r1,'b.');
plot(sgm4_lon_r2,sgm4_lat_r2,'ro');
title('Segment 1-row1  and Segment 4-row 2'); 

subplot(2,2,3);
hold on;
plot(sgm2_lon_r2,sgm2_lat_r2,'b.');
plot(sgm3_lon_r1,sgm3_lat_r1,'ro');
title('Segment 2-row2  and Segment 3-row 1'); 

subplot(2,2,4);
hold on;
plot(sgm2_lon_r1,sgm2_lat_r1,'b.');
plot(sgm3_lon_r2,sgm3_lat_r2,'ro');
title('Segment 2-row1  and Segment 3-row 2'); 

% Cut topography from GLBb into ARCc grid:
Hg2a=HH*nan;
% Part that does not need any rotation/swapping:
dmm=HG(j1:end,i1:i2);
[a1,a2]=size(dmm);
Hg2a(1:a1,1:a2)=dmm;
[a1,a2]=size(dmm);
%ipv=i3+nn/2-1;  
% Segment3 -> Segment 2
dmm=HG(j3:end-2,i3:ipv);
dmm=fliplr(dmm);
dmm=flipud(dmm);
Hg2a(a1+1:end,nn/2+1:nn)=dmm;

% Segment4 -> Segment 1:
dmm=HG(j3:end-2,ipv+1:i4);
dmm=fliplr(dmm);
dmm=flipud(dmm);
Hg2a(a1+1:end,1:nn/2)=dmm;

I=find(isnan(Hg2a));
if ~isempty(I);
  fprintf('GLB -> ARCc topo error \n');
  error(' Not all topo points have been feeled')
end

% Topo difference:
D=Hg2a-HH;

figure(1); clf;
axes('position',[0.05 0.05 0.4 0.9]);
HH(HH>0)=nan;
pcolor(HH); shading flat;
axis('equal');
set(gca,'xlim',[1 nn],'ylim',[1 mm]);
caxis([-6000 0]);
title('ARCc topo');

axes('position',[0.52 0.05 0.4 0.9]);
Hg2a(Hg2a>0)=nan;
pcolor(Hg2a); shading flat;
axis('equal');
set(gca,'xlim',[1 nn],'ylim',[1 mm]);
caxis([-6000 0]);
title('GLBb topo in ARC');



%txtbtm='/hycom_arc08/NOPP_rivers/compare_T07_T09.m';
txtbtm='/hycom_arc08/NOPP_rivers/topo_GLBbT07_ARCcT09.m';
bottom_text(txtbtm);













