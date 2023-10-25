% Correct ARCc0.04 topo 17
% Open up the strait between 
% the Gulf of Boothia and 
% Fox Basin and make Baffin Island as Island
% not peninsula
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close

f_nc=1; % save also netCDF bath. file;
f_cice_grid = 0; % =1 - prepare cice grid that matches corrected TOPO
                 % suggested: use make_cicegrid.csh that calls Alan's
		 % routines, as my code has not been validated 

RG='ARCc0.04'; % 
hg=2^100;

PTH.topo=sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/%s/topo_grid/',RG);
T = 17;     % topo #
TV='17DD';  % new name 
pthT = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/%s/topo_grid/',RG);
pthout=pthT;



% Topo/grid files, input:
flda=sprintf('%sdepth_%s_%2.2i.a',pthT,RG,T);
%flda=sprintf('%sdepth_ARCc0.04_11D.a',pthT);
flga=sprintf('%sregional.grid.a',pthT);
flgb=sprintf('%sregional.grid.b',pthT);

% Read in topo from *a file:
f1=fopen(flgb);  % read I,J from regional.grid.b
aa=fgetl(f1);
dmm=aa(2:8);
IDM=str2num(dmm);
aa=fgetl(f1);
dmm=aa(2:8);
JDM=str2num(dmm);
IJDM=IDM*JDM;
fclose(f1);

npad=4096-mod(IJDM,4096);

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % read lon/lat from regional grid file
grid_fid=fopen(flga,'r');

[plon,count]=fread(grid_fid,IJDM,'float32','ieee-be');
stat=fseek(grid_fid,4*(npad+IJDM),-1);
if stat<0
  error('Reading grid file ...');
end
[plat,count]=fread(grid_fid,IJDM,'float32','ieee-be');

disp('Reading lat/lon  ...')
plon=(reshape(plon,IDM,JDM))';
plat=(reshape(plat,IDM,JDM))';

fclose(grid_fid);

% --------------------------
I=find(plon>180);
plon(I)=plon(I)-360;
I=find(plon<-180);
plon(I)=plon(I)+360;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read bathymetry from regional.depth.a

depth_fid=fopen(flda,'r');

%fseek(depth_fid,6*4*(npad+IJDM),-1) % this is confusing, should not be needed
[h,count]=fread(depth_fid,IJDM,'float32','ieee-be');
%y=find(h>1e10);
%h(y)=nan;
HH=reshape(h,IDM,JDM)';
clear h;
fclose(depth_fid);

HH(HH>1e10)=nan;
HH=-HH;
HH(HH>-0.1)=0;
HH(isnan(HH))=100;
[m,n] = size(HH);


% Fix:
% make sure the end points are in the ocean
i11=759;
i22=766;
j11=1882;
j22=1897;
HH11=HH;
for ii=i11:i22
  h=HH(j11:j22,ii);
  jland1 = min(find(h>=0))-1;
  jland2 = max(find(h>=0))+1;
  npj=jland2-jland1;
  dh=(h(jland2)-h(jland1))/npj;
  h11=h(jland1)+dh*(1:npj-1)';
  h(jland1+1:jland2-1)=h11;
  HH11(j11:j22,ii)=h;
end


% ========================
%
%   Write new topo to ARC"c"
%
% =======================
fina=sprintf('%sdepth_%s_%s.a',PTH.topo,RG,TV);
finb=sprintf('%sdepth_%s_%s.b',PTH.topo,RG,TV);
fida=fopen(fina,'w','ieee-be');
fidb=fopen(finb,'wt');
fprintf('Writing new topo %s\n',fina);

% Write 5-line header
fprintf(fidb,'bathymtry 30sec GEBCO_08 20091120 data, 5x5 avrg; Fixed: closed strait G. of Boothia\n');
fprintf(fidb,'i/jdm = %i %i; plon,plat range = %8.5f %8.5f %8.5f %8.5f\n',...
	n,m,min(min(plon)),max(max(plon)),min(min(plat)),max(max(plat)));
fprintf(fidb,'Filled single-width inlets and (B&C grid) enclosed seas except the Black Sea.\n');
fprintf(fidb,'Coastline at 0.1m.  Smoothed 1x.  S.Ocean under Antarctic ice shelves\n');
fprintf(fidb,'subregion of GLBc0.04-17; landmask all edges, Svalbard, Med. and Black Sea\n');
fprintf(fidb,'min,max depth = %11.5f %11.5f\n',...
	abs(max(max(HH11(HH11<0)))),abs(min(min(HH11))));

fclose(fidb);
HH11(isnan(HH11))=100;

% Land -> huge #
% all depths >0
Ioc=find(HH11<0);
Ilnd=find(HH11>=0);
HH11=abs(HH11);
HH11(Ilnd)=hg;

IDM=n;
JDM=m;
IJDM=IDM*JDM;
npad=4096-mod(IJDM,4096);
toto=ones(npad,1);
A=reshape(HH11',IJDM,1);
fwrite(fida,A,'float32');
fwrite(fida,toto,'float32');
fclose(fida);

fprintf('New bathymetr is done %s\n',fina);

if f_cice_grid > 0;

% ------------------------
% 
% Correct sea ice topo file:
% make_cicegrid.csh  
%
% Use script and Alan's code
% or run this piece of matlab code:
%
% read bathymetry from cice.regional.r
%
% in grid2cice.f (located at:
% /pacific/abozec/HYCOM/hycom/ALL2/cice/src
% ice grid is written as direct access array:
% each record is IDMo x JDMo double precision (real*8)
%c        kmt    land mask array (0,1)
%c        ulati  latitude  of u-cell centers (radians)
%c        uloni  longitude of u-cell centers (radians)
%c        htn    length of northern edge of t-cell (m)
%c        hte    length of eastern  edge of t-cell (m)
%c        anglet conversion on t-cell between cice and lat-long grids (radians)
%c        tlati  latitude  of t-cell centers (radians)
%c        tloni  longitude of t-cell centers (radians)
% 
% ------------------------
  flrga = sprintf('%sregional.cice_%2.2i.a',pthT,T);
  flrgb = sprintf('%sregional.cice_%2.2i.b',pthT,T);
  flrgr = sprintf('%sregional.cice_%2.2i.r',pthT,T); % old cice.r
%  flrgr = sprintf('%sregional.cice.T%s.r',pthT,TV);% New created by grid2cice.f Alan
  flrgaN = sprintf('%sregional.cice_%s.a',pthT,TV); % new created here
  %flrgbN = sprintf('%sregional.cice_%s.b',pthT,TV); % nothing will change here
  flrgrN = sprintf('%sregional.cice_%s.r',pthT,TV);

  Lmsk = HH*0;
  Lmsk(isnan(HH11))=0;
  Lmsk(HH11<0)=1;

  % Fix regional.cice*.r
  depth_fido=fopen(flrgr,'r');
  fidr_out = fopen(flrgrN,'w');

  dmm=fread(depth_fido,1,'float64','ieee-be');
  kmt=fread(depth_fido,[IDM,JDM],'float64','ieee-be');
  kmt=kmt';
  dmm=fread(depth_fido,1,'float64','ieee-be');

  % Need to shift left by 1/2 grid due to different staggerd grids
  % used in CICE and HYCOM
  kmtN = kmt;
  kmtN(:,2:end-1)=Lmsk(:,3:end);
  %dL = Lmsk-kmtN;
  %dC = kmt-kmtN; % changes in the cice mask

  kmtN = kmtN';
  
  fwrite(fidr_out,dmm,'float64');
  fwrite(fidr_out,kmtN,'float64');
  fwrite(fidr_out,dmm,'float64');

  error('The pice of code is not finished, need to read/write other fields...');

  fclose(depth_fido);
  fclose(fidr_out);

end;



if f_nc==0, return; end;

flda=fina;
flga=sprintf('%sregional.grid.%s.a',PTH.topo,RG);
flgb=sprintf('%sregional.grid.%s.b',PTH.topo,RG);
pthout=PTH.topo;
foutB=sprintf('depth_%s_%s',RG,TV);
js=1;
je=m;
is=1;
ie=n;
ref1 = 'NRL SSC, bath. from 30-scnd GEBCO_08';
ref2 = 'edited COAPS FSU D.Dukhovskoy';
sub_convTOPO2nc(RG,TV,flda,flga,flgb,pthout,foutB,js,je,is,ie,ref1,ref2);



