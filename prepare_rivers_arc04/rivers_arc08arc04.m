% Rivers for Topography T17DD to ARCc0.04
% Using ARCc0.08 T11D
% Units in hycom river*.[ab] is  m/s
%c --- initialize input of river (precip bogas) forcing field
%c --- units of rivers are m/s    (positive into ocean)
%c --- rivers is always on the p grid.
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
startup;

close all
clear

R = 'ARCc0.04';
E = '010';

AA.Region1   = 'ARCc0.08';
AA.Region2   = 'ARCc0.04';
AA.Expt1     = '110';
AA.Expt2     = E;
AA.Topo1     = 11;
AA.Topo2     = 17;
AA.pthdat1   = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/force/rivers/';
AA.pthtopo1  = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
AA.pthdat2   = '/Net/mars/ddmitry/hycom/ARCc0.04/force/rivers/';
AA.pthtopo2  = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/topo_grid/';

%YY=2004;
% ARCc0.08 rivers, 
flriv1a = sprintf('%srivers_%2.2i.a',AA.pthdat1,AA.Topo1);
flriv1b = sprintf('%srivers_%2.2i.b',AA.pthdat1,AA.Topo1);
fltopo1 = sprintf('%sdepth_ARCc0.08_%2.2i.nc',AA.pthtopo1,AA.Topo1);

flriv2a = sprintf('%srivers_%2.2iDD.a',AA.pthdat2,AA.Topo2); 
flriv2b = sprintf('%srivers_%2.2iDD.b',AA.pthdat2,AA.Topo2);
fltopo2 = sprintf('%sdepth_ARCc0.04_17DD.nc',AA.pthtopo2); % corrected topo - open strait

% Old topo:
HH1  = nc_varget(fltopo1,'Bathymetry');
alat = nc_varget(fltopo1,'Latitude');
elon = nc_varget(fltopo1,'Longitude');
LAT1  = alat;
LON1  = elon;
[mm1,nn1]= size(HH1);
IDM1=nn1;
JDM1=mm1;
IJDM1=IDM1*JDM1;
npad1=4096-mod(IJDM1,4096);
toto1=ones(npad1,1);

[DX1,DY1]=sub_dx_dy(LON1,LAT1);
ACell1=DX1.*DY1;

% Get new topo and grid:
HH2  = nc_varget(fltopo2,'Bathymetry');
LAT2 = nc_varget(fltopo2,'Latitude');
LON2 = nc_varget(fltopo2,'Longitude');
[mm,nn]= size(HH2);
[m,n]= size(HH2);
IDM=n;
JDM=m;
IJDM=IDM*JDM;
npad=4096-mod(IJDM,4096);
toto=ones(npad,1);
toto2 = toto;

[DX2,DY2]=sub_dx_dy(LON2,LAT2);
ACell2=DX2.*DY2;

%keyboard

% ==============================
%
%      CORRECT/RELOCATE RIVERS and 
%      WRITE OUTPUT FILES
%
% ==============================

% Open files with Rivers for IO:
faold = fopen(flriv1a,'r','ieee-be');
fbold = fopen(flriv1b,'r');
fanew = fopen(flriv2a,'w');
fbnew = fopen(flriv2b,'wt');


fprintf('Opening files to read: \n');
fprintf('%s\n',flriv1a);
fprintf('%s\n',flriv1b);
fprintf('Opening files to write: \n');
fprintf('%s\n',flriv2a);
fprintf('%s\n',flriv2b);

% Write heading:
for nl=1:5
  aa=fgetl(fbold);
  disp(aa);

  if nl==4
    aa=sprintf('Rivers ARCc0.08 interpoltd to ARCc0.04 T17 correctd CAA strght');
  end

  fprintf(fbnew,[aa,'\n']);
end


for k=1:12
  fprintf('Month = %i\n',k);
% Read *.b:  
  aa=fgetl(fbold);
  disp(aa);

% Read HYCOM runoff:
  R1=fread(faold,IJDM1,'float32'); % read 2D field
  dm1=fread(faold,npad1,'float32');  % Padding = size(toto)
  if size(dm1) ~= size(toto1)
    error('Padding in HYCOM file ???');
  end
%  toto=dm1;
  clear dm1

%  I=find(A>1e10);  % there are no land masks in rivers
%  A(I)=NaN;
  R1 = reshape(R1,IDM1,JDM1)';
  R2 = sub_riversARCc08toARCc04(HH1,R1,ACell1,LON1,LAT1,...
				HH2,ACell2,LON2,LAT2);
  
  
  

  IR=find(R2>0 & HH2>=0);
  if ~isempty(IR),
    error('Found river grid cells on land\n');
  end
%keyboard  

  
  
% ----------------
  f_chck=0;
  if f_chck==1
    figure(10); clf;
    contour(HH2,[0.1 0.1],'k','Color',[0 0.8 0.8]);
    hold on;
%    contour(HHo,[0.1 0.1],'k','Color',[0.8 0.8 0.8]);
    R=R2;
    R(R==0)=nan;
    pcolor(R); shading flat;
    caxis([0 5e-5]);
    keyboard  
  end
% ----------------

% Check river points

  fprintf('Writing HYCOM files %s\n',flriv2a);  
  fprintf('Writing HYCOM files %s\n\n',flriv2b);  

  dmm=R2';
  dmm=reshape(dmm,IJDM,1);
  minh=min(dmm);
  maxh=max(dmm);
  is=strfind(aa,'=');
  ch=aa(1:is);
  anew=sprintf('%s%3i    %10.7E   %10.7E',ch,k,minh,maxh);
  fprintf(fbnew,[anew,'\n']);
  fwrite(fanew,dmm,'float32','ieee-be');
  fwrite(fanew,toto,'float32','ieee-be');
%keyboard
end;  % for k - months

fclose(fanew);
fclose(fbnew);
fclose(faold);
fclose(fbold);
















