% Convert Rivers from Topography T09 to T11
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

PTH.data='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/force/rivers/';
PTH.topo='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
%flt=[pth,'regional.grid.b'];

ntopo1=09;
ntopo2=11;

%YY=2004;
% Old river runoff files: use corrected runoff
flriva=sprintf('%srivers_%2.2i_corrected.a',PTH.data,ntopo1);
flrivb=sprintf('%srivers_%2.2i_corrected.b',PTH.data,ntopo1);
%flriva=sprintf('%srivers_09_Greenland_%i.a',PTH.data,YY);
%flrivb=sprintf('%srivers_09_Greenland_%i.b',PTH.data,YY);
fltopo_old=sprintf('%sdepth_ARCc0.08_%2.2i.nc',PTH.topo,ntopo1);

fltopo_new=sprintf('%sdepth_ARCc0.08_%2.2i.nc',PTH.topo,ntopo2);
flriv_new_a=sprintf('%srivers_%2.2i.a',PTH.data,ntopo2); 
flriv_new_b=sprintf('%srivers_%2.2i.b',PTH.data,ntopo2);

% Old topo:
HHo  = nc_varget(fltopo_old,'Bathymetry');

% Get new topo and grid:
HH   = nc_varget(fltopo_new,'Bathymetry');
alat = nc_varget(fltopo_new,'Latitude');
elon = nc_varget(fltopo_new,'Longitude');
LAT  = alat;
LON  = elon;
[mm,nn]= size(HH);
[m,n]= size(HH);
IDM=n;
JDM=m;
IJDM=IDM*JDM;
npad=4096-mod(IJDM,4096);
toto=ones(npad,1);

[DX,DY]=sub_dx_dy(LON,LAT);
ACell=DX.*DY;


% ==============================
%
%      CORRECT/RELOCATE RIVERS and 
%      WRITE OUTPUT FILES
%
% ==============================

% Open files with Rivers for IO:
faold = fopen(flriva,'r','ieee-be');
fbold = fopen(flrivb,'r');
fanew = fopen(flriv_new_a,'w');
fbnew = fopen(flriv_new_b,'wt');

% Write heading:
for nl=1:5
  aa=fgetl(fbold);
  disp(aa);

  if nl==4
    aa=sprintf('Rivers corrected from T09 relocated to T11');
  end

  fprintf(fbnew,[aa,'\n']);
end


for k=1:12
% Read *.b:  
  aa=fgetl(fbold);
  disp(aa);

% Read HYCOM runoff:
  A=fread(faold,IJDM,'float32'); % read 2D field
  dm1=fread(faold,npad,'float32');  % Padding = size(toto)
  if size(dm1) ~= size(toto)
    error('Padding in HYCOM file ???');
  end
%  toto=dm1;
  clear dm1

%  I=find(A>1e10);  % there are no land masks in rivers
%  A(I)=NaN;
  A=reshape(A,IDM,JDM)';
  
  R=A;
  A = sub_rivers2coast(HH,R);
  IR=find(A>0 & HH>=0);
  if ~isempty(IR),
    error('Found river grid cells on land\n');
  end
%keyboard  

  
  
% ----------------
  f_chck=0;
  if f_chck==1
    figure(10); clf;
    contour(HH,[0.1 0.1],'k','Color',[0 0.8 0.8]);
    hold on;
    contour(HHo,[0.1 0.1],'k','Color',[0.8 0.8 0.8]);
    R=A;
    R(R==0)=nan;
    pcolor(R); shading flat;
    caxis([0 5e-5]);
    keyboard  
  end
% ----------------

% Check river points

  fprintf('Writing HYCOM files %s\n',flriv_new_a);  
  fprintf('Writing HYCOM files %s\n\n',flriv_new_b);  

  dmm=A';
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











