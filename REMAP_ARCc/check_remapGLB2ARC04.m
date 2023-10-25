% Check created archv files on ARCc grid
% from GLBb - 1st step, no interpolation yet
% simply subsmaple GLBb and rotate to ARCc:
%
% Original files: GLBc0.04 GOFS3.5 analysis
%
% In the ARCc - the half that is 
% rotated from the GLBb grid (up from ARCc j = xxx)
% all vectors (U and V components) need to be rotated by 180degrees
% to match + direction in the ARCc grid
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close

R    = 'ARCc0.04';
%fld  = 'salin'; % field to plot
fld  = 'u-vel.'; % field to plot
%fld  = 'v_btrop'; % field to plot
lplt = 1;       % layer to plot
YR   = 2017;
DD   = 2;
EE   = 216; % new reanalysis
expt = '21.6';

fprintf('Field: %s\n',fld);


switch(fld),
 case('salin')
  c1 = 32; 
  c2 = 34;
 case ({'u-vel.','v-vel.'})
  c1 = -0.4;
  c2 = 0.4;
 case ({'u_btrop','v_btrop'})
  c1 = -0.2;
  c2 = 0.2;
end


pathglb = '/Net/kronos/ddmitry/hycom/GOFS3.5/';
patharc = '/Net/kronos/ddmitry/hycom/GOFS3.5/';
pathtopo= '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/topo/';
pathgmap = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/topo/';
fnmGMAP  = 'regional.gmapi_GLBc0.04';
fnmGLB   = sprintf('%3.3i_archv.%i_%3.3i_00',EE,YR,DD);
fnmARC   = sprintf('archv_arcT17L41.%i_%3.3i_00_stp1',YR,DD);  


%IDMg     = 4500;
%JDMg     = 3298;
%IJDMg = IDMg*JDMg;
IDMa     = 3200;
JDMa     = 5040;
IJDMa    = IDMa*JDMa;

% indices of GLB that correspnd ARCc
iArc1=2449;
iArc2=4048;
iArc3=453;
iArc4=2052;
jArc1=2051;
jArc2=2051;
jArc3=2025;
jArc4=2025;

pthtglb = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/topo_grid/';
pthtarc = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/topo_grid/';

hg=2^100;
rg=9806;  % convert pressure to depth, m

yr   = 2001;
iday = 36;
hr   = 0;
TV   = '07'; % arc08 & global

% Read global fields
% Note that for new reanalysis global Topo = 11
% In GLBb archm u-vel and v-vel are total velocities
% where as in nest files archv - these are baroclinic anomalies
frga   = sprintf('%sregional.grid_GLBb0.08_07.a',pthtglb);
frgb   = sprintf('%sregional.grid_GLBb0.08_07.b',pthtglb);
fdpthGa = sprintf('%sdepth_GLBb0.08_07.a',pthtglb);
fdpthGb = sprintf('%sdepth_GLBb0.08_07.b',pthtglb);

fidRGb = fopen(frgb,'r');  % read I,J from regional.grid.b
aa  = fgetl(fidRGb);
dmm = aa(2:8);
IDMg = str2num(dmm);
aa = fgetl(fidRGb);
dmm = aa(2:8);
JDMg = str2num(dmm);
IJDMg = IDMg*JDMg;
fclose(fidRGb);
npadg=4096-mod(IJDMg,4096);
fprintf('Global: IDM=%i, JDM=%i\n',IDMg,JDMg);

fprintf('Reading Global topo %s\n',fdpthGa);
% read lon/lat from GLBb regional grid file
fidRGa = fopen(frga,'r');
[plon,count] = fread(fidRGa,IJDMg,'float32','ieee-be');
fseek(fidRGa,4*(npadg+IJDMg),-1);
[plat,count] = fread(fidRGa,IJDMg,'float32','ieee-be');

fidtGa = fopen(fdpthGa,'r');
dmm = fread(fidtGa,IJDMg,'float32','ieee-be');
HHg = reshape(dmm,[IDMg,JDMg])';
I = find(HHg>1e10);
HHg = -HHg;
HHg(I)= 100;

%figure(1); clf;
%contour(HHg,[0 0],'k');

% Plot fields:

finga = sprintf('%s%s.a',pathglb,fnmGLB);
fingb = sprintf('%s%s.b',pathglb,fnmGLB);
[F,n,m,l] = read_hycom(finga,fingb,fld,'r_layer',lplt);

Ag = squeeze(F(1,jArc3:end,iArc3:iArc2));
clear F
Ag(Ag>1e10)=nan;
Hg=HHg(jArc3:end,iArc3:iArc2);

figure(1); clf;
pcolor(Ag); shading flat;
hold on;
contour(Hg,[0 0],'k');
caxis([c1 c2]);
colorbar;
stt=sprintf('%s, Layer: %i, GLBb0.08, T07, %s',fld,lplt, fnmGLB);
sprintf(stt);
axis('equal');
title(stt,'Fontsize',14,'Interpreter','none');




fprintf('Reading ARCc %s\n',fnmARC);
finaa = sprintf('%s%s.a',pathtopo,fnmARC);
finab = sprintf('%s%s.b',pathtopo,fnmARC);
[F,n,m,l] = read_hycom(finaa,finab,fld,'r_layer',lplt);

Aa = squeeze(F(1,:,:));
Aa(Aa>1e10)=nan;

figure(2); clf;
pcolor(Aa); shading flat;
hold on;
contour(HHa,[0 0],'k');
caxis([c1 c2]);
colorbar;
stt=sprintf('%s, Layer: %i, ARCc0.04, T27, %s',fld,lplt,fnmGLB);
sprintf(stt);
axis('equal');
set(gca,'xlim',[1 IDMa],'ylim',[1 JDMa]);
title(stt,'Fontsize',14,'Interpreter','none');


%
% Read ARCc0.04:
%fltopa=sprintf('%sdepth_ARCc0.04_17.a',pathtopo);
%fltopb=sprintf('%sdepth_ARCc0.04_17.b',pathtopo);
fltopa=sprintf('%sdepth_ARCc0.04_GOFS3.5_27.a',pathtopo);
fltopb=sprintf('%sdepth_ARCc0.04_GOFS3.5_27.b',pathtopo);
flgrd = sprintf('%sregional.grid.ARCc0.04',pathtopo);

GRD = read_grid_bath(flgrd,fltopa);
HH = GRD.Topo;
I=find(HH>1e20);
HH=-HH;
HH(I)=100;

% Check matching land points:
Il1=find(HH>=0);
Il2=find(isnan(Aa));

if (length(Il1)==length(Il2))
  fprintf('%s matches ARCc0.04\n',fltopa);
else
  fprintf('%s DOES NOT MATCH  ARCc0.04\n',fltopa);
end 

% Plot T27 from GLBc0.04 --> ARCc0.04 grid
%figure(2); clf;



