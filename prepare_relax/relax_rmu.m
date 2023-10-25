% For relaxation or nesting, need to specify relaxation weights
% within the buffer zone:
%  relax_rmu.[ab] - for relaxation
%                   this is needed if there is 3D relaxation
%                   to climatology
%  nest_rmu.[ab] - nesting, can be the same as relax_rmu.[ab]
%                  needed if model is nested in the outer model
%                  these are relaxation weights within the 
%                  buffer zone
%
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close

ptharc  = '/Net/kronos/ddmitry/hycom/ARCc0.08/tmp_files/';
%pthglb  = '/Net/kronos/ddmitry/hycom/GLBb0.08/expt_69.1/restart/';
pthin  = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/relax/010/';
pthout = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/relax/110/41layers_T11/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';

ntopo1=09;
ntopo2=11;

% Read in sigma layers:
fina = sprintf('%srelax_rmu.a',pthin);
finb = sprintf('%srelax_rmu.b',pthin);
fouta = sprintf('%srelax_rmu_T%2.2i.a',pthout,ntopo2);
foutb = sprintf('%srelax_rmu_T%2.2i.b',pthout,ntopo2);

fida = fopen(fina,'r','ieee-be');
fidb = fopen(finb,'r');

aa='a';
pat='[ij]dm';
while isstr(aa)
  aa = fgetl(fidb);
  if ~isstr(aa); break; end;
  ii=regexp(aa,pat);
  if ~isempty(ii),
    ind=findstr(aa,'=');
    ch=aa(ind+1:end);
    dmm=sscanf(ch,'%i');
    IDM=dmm(1);
    JDM=dmm(2);
    break
  end
end
fclose(fidb);
      
IJDM=IDM*JDM;
npad=4096-mod(IJDM,4096);
toto=ones(npad,1);
kold=32;
knew=41;

fprintf('ARCc domain, ID=%i JD=%i\n',IDM,JDM);

% Read relaxaion file:
fida = fopen(fina,'r');
A = fread(fida,IJDM,'float32','ieee-be');
A = (reshape(A,IDM,JDM))';
fclose(fida);


% Get new topo and grid:
fltopo_new=sprintf('%sdepth_ARCc0.08_%2.2i.nc',pthtopo,ntopo2);
HH   = nc_varget(fltopo_new,'Bathymetry');
alat = nc_varget(fltopo_new,'Latitude');
elon = nc_varget(fltopo_new,'Longitude');
LAT  = alat;
LON  = elon;
[mm,nn]= size(HH);


fltopo_old=sprintf('%sdepth_ARCc0.08_%2.2i.nc',pthtopo,ntopo1);
HHo  = nc_varget(fltopo_old,'Bathymetry');

f_chngR=0; % change relaxation
if f_chngR>0
  nrx=20;
  rmn=0.04;  % strongest relaxation
  rmx=10;   % weakest relaxation
  xx=(1:nrx);
  a2=rmn;
  a1=(rmx-a2)/(nrx-1)^2;
  rlx_day=a1*(xx-1).^2+a2;  % qudratic decay of relaxation time in the buffer zone
  %rlx_day=(rmn:(rmx-rmn)/(nrx-1):rmx);  % relaxation of days, e-folding time; lines decay
  rlx=1./(rlx_day*86400);
  plot(rlx,'.-');
  xlabel('Grid points in buffer zone');
end;

RMU = zeros(JDM,IDM);
% Copy relaxation from T09 file:
npp=50;
% Southern OB (Atlantic):
rlx=A(1:npp,400);
nrx=length(find(rlx>0));
for j=1:npp
  RMU(j,:)=rlx(j);
end

% NOrhtern(Pacific)
rlx=A(:,400);
for j=JDM-npp:JDM
  RMU(j,:)=rlx(j);
end

% Exclude land
I=find(HH>0);
RMU(I)=0;

f_chck=0;
if f_chck>0
  cmp=colormap('jet');
  cmp(1,:)=[1 1 1];

  figure(1); clf;
  pcolor(RMU); shading flat;
  colormap(cmp);
  hold
  contour(HH,[0 0],'r');
  axis('equal');
  set(gca,'xlim',[1 IDM],...
	  'ylim',[1 JDM]);
  title('Relaxation Zone');
end;

rr1=min(min(RMU(RMU>0)));
rr2=max(max(RMU));
rmn=(1/rr2)*1/86400; % days, min relax
rmx=(1/rr1)*1/86400; % days, max relax

% Write relax files:
fid_aout = fopen(fouta,'w');
fid_bout = fopen(foutb,'wt');

fprintf(fid_bout,'Relaxation Mask, ARCc0.08 Topo 11\n');
fprintf(fid_bout,...
 'S(Atl) and N(Pacif) boundaries: %i gridpts with %5.3f-%5.3f day e-folding time \n',nrx,rmn,rmx);
fprintf(fid_bout,' \n\n');
fprintf(fid_bout,'i/jdm = %i  %i\n',IDM,JDM);
fprintf(fid_bout,'     rmu: range =      %d  %d\n',...
      min(min(RMU)), max(max(RMU)) );
fclose(fid_bout);

RMU=RMU';
A=reshape(RMU,IJDM,1);
fwrite(fid_aout,A,'float32','ieee-be');
fwrite(fid_aout,toto,'float32','ieee-be');  % padding at the end
fclose(fid_aout);

fprintf('Files are in: %s \n',pthout);
