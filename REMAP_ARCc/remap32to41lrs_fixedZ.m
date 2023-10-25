% Test code for remapping
% 32 v layers from GLBb0.08 T07 experiment 19.0
% onto 41 v layers ARCc0.07 T11
% assuming top new layers are fixed-depth layers
%
% The first 2 steps (GLB-> ARC and ARC T07 ->T11) 
% should be already done
% This code uses ARCc T11 32 layer archive files
% created from GLBb0.08 T07
% and interpolates into 41 layers used in GOFS3.1
%
% Use layer depths from PHC climatology file
% for 41 layers (see: RELAX_PHC/*)
%
% Targ. densities of these 41 layers are such that
% added layers are at the top of the ocean
% above sigma2=35.50 (L41=24) all layers beneath it
% match 32-layer grid
% The added layers in the upper 100-200 m
% are anticipated to be FIXED depth layers 
%
% Inspection of hybrid layers generated for HPC climatology
% indicates that the 2 last layers above sigma2=35.50 may be 
% not fixed-depth
% 
% This algorithm assumes:
% 1) All layers L41<24 (sigma2=35.50) are fixed Z/sigma depths
% 2) Depth of the bottom interface of L41(24) = L32(14) 
%    i.e. sum[i=1,24]L41(i)=sum[i=1,14]L32(i), however 
%    at some locations 32-layer 14 is shallower than
%    L41(24) and even L41(23) - will need to adjust this
% 
% 3) The bottom layer L32(32), sigma2=37.48 does not exist in L41
%    will need to merge the last two L32 layers in 1 for L41

addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close

f_intrp=1; % 1- interpolate; 
           % 2 - load saved mat with interp. flds, write *a,*b

pthbin  = '/Net/kronos/ddmitry/hycom/ARCc0.08/tmp_files/';
pthout  = '/Net/kronos/ddmitry/hycom/ARCc0.08/archv_41layers/';
pthrlx = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/relax/110/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';

hg=2^100;
rg=9806;  % convert pressure to depth, m

yr=1993;
iday=1;
hr=0;
TV=11;
%nsigma=14; % # of sigma levels from blkdat.input
%mo=7; % month, 1 or 7

% Topo:
fltopo=sprintf('%sdepth_ARCc0.08_%2.2i.nc',pthtopo,TV);
HH  = nc_varget(fltopo,'Bathymetry');
LAT = nc_varget(fltopo,'Latitude');
LON = nc_varget(fltopo,'Longitude');
[mm,nn] = size(HH);
%HH(HH>0)=nan;


IDM=nn;
JDM=mm;

% Array of z/s depths for not sigma region:
% Read dp0k from blkdat.input
% with specified z-level spacing
pthb='/home/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers/';
fnm=sprintf('%sblkdat.input_ARCc0.08_41lev',pthb);
VL=read_layers_blkdat(fnm);
dPdeep=VL.dp0k;
nsigma = VL.Nsigma_layers; % sigma layers on shelf

% Read dP41 layer thicknesses from
% climatology:
% This fields are needed for remapping
% for searching sigma-layer regions (could
% be done based on blkdat.input dp0k spacing)
% and for checking local       bottom depths
% and Lx_new layer depth
mo=1;
finRa = sprintf('%s41layers_T%2.2i/relax.m%2.2i.a',pthrlx,TV,mo);
finRb = sprintf('%s41layers_T%2.2i/relax.m%2.2i.b',pthrlx,TV,mo);
fld='thknss';
[Fg,nr,mr,lr] = read_hycom(finRa,finRb,fld);
lr_new=size(Fg,1);

Fg(Fg>1e10)=nan;
Fg(Fg<0.1)=nan;

Fg=Fg./rg;
DP_new=Fg;

Lx_new=24; % Layer41 where 41 and 32-lr grids should match
Lx_old=14; % layer from 32-layer where 2 grid match
%keyboard

% Interf. depths and depths of the middle of the layers (m)
% NOTE: sign convention: depths are negative
% dP - Pa or m
[ZZn1,ZMn1] = sub_thck2dpth(DP_new); % 


% Find points not in sigma-region, i.e.
% where depth> sum over (nsigma dP)
sgmDP=squeeze(nansum(DP_new(1:nsigma,:,:),1));
sgmDP(sgmDP==0)=nan;
HH(HH>0)=nan;
sgmH  = abs(1-abs(sgmDP)./abs(HH)); % sigma-coord. points
Isgm=find(sgmH<1e-3); % sigma points
%[js,is]=ind2sub(size(HH),Isgm);
sgmH(Isgm)=nan;
Insgm = find(~isnan(sgmH)); % not-sigma coord. points

dP=squeeze(DP_new(Lx_new-1,:,:));
dP(Isgm)=nan;

% Plot the next to fixed layer depth (where 32 and 41 grids match)
% - should be
% constant throughout the domain
% It is not, right at the boundaries
% assume the depth is constant
HH(isnan(HH))=100;
f_pltdp=0;
if f_pltdp>0
  figure(1); clf;
  contour(HH,[0 0],'k');
  hold on;
  pcolor(dP); shading flat;
end;


% for checking
% Use some Arctic Ocean location
%iZ=800;
%jZ=1500;
%dPdeep=squeeze(DP_new(1:Lx_new-1,jZ,iZ));
Zdeep=-cumsum(dPdeep);

% On sigma-regions (shelves) - sigma coordinates
% Below Lx_new (41 layers) - layers match, except for the bottom layer

% Interpolate ARCc T11 32 layers -> ARCc T11 41 layers:
% Because of space limitation, follow the algorithm:
% interpolate by variables and save in *.mat files
% For output, go by layers and read in variable, layer
% write the output

% Input:
finaOLD = sprintf('%sarchv_arcT%2.2i.%4.4i_%3.3i_%2.2i.a',pthbin,TV,yr,iday,hr);
finbOLD = sprintf('%sarchv_arcT%2.2i.%4.4i_%3.3i_%2.2i.b',pthbin,TV,yr,iday,hr);
fld='thknss';
[Fg,nr,mr,lr] = read_hycom(finaOLD,finbOLD,fld);
lr_old=size(Fg,1);

Fg(Fg>1e10)=nan;
Fg(Fg<0.1)=nan;

Fg=Fg./rg;
DP_old=Fg;
[ZZo,ZMo] = sub_thck2dpth(DP_old); % 
clear Fg

% First: update new vertical grid
% above Layer 24 - fixed depths
% below - identically matches the 32-layer grid
% the bottom-most layer combines 2 bottom layers of 32 grid
%keyboard
if f_intrp==1
  ZZn=sub_adjust_vlayers(ZZn1,ZZo,Lx_new,Lx_old,HH,nsigma,dPdeep);
  fmat=sprintf('%sZZ_interp32to41.mat',pthout);
  fprintf('Saving adjusted v.layers %s\n',fmat);
  save(fmat,'ZZn');
else
  fmat=sprintf('%sZZ_interp32to41.mat',pthout);
  load(fmat);
end  
[nzN,mm,nn]=size(ZZn);
nlN=nzN-1;
clear ZZn1 ZMn1 ZMn ZMo 


f_chck=0; % check grid
if f_chck>0
  vname='temp';
  [Fold,n,m,l] = read_hycom(finaOLD,finbOLD,vname);% read in 32-layer field
  Fold(Fold>0.0001*hg)=nan;
  Fnew=zeros(nlN,mm,nn);
%  IJ=[200,15;1020,15];
  IJ=[1,2518;1600,2518];
  sub_plot32v41(Fnew,Fold,HH,IJ,ZZn,ZZo,vname);
end

% If memory not an issue
% interpolate all 3D fields then 
% write into *a file
% otherwise will have to do 1 variable at a time
VARS={'temp','salin','u-vel.','v-vel.'};
nv=length(VARS);

%for testing - leave only relaxation zones:
%HHm=HH;
%HHm(30:mm-30,:)=100;

if f_intrp==1
for kk=1:nv;
  tic;
%vname = 'temp';
  vname=VARS{kk};
  Fnew = sub_interp_fld32to41(vname,finaOLD,finbOLD,...
			      HH,ZZn,ZZo,nsigma,Lx_new, Lx_old);
  fprintf('%s Interpolation done ...\n\n',vname);
  switch(vname);
   case('temp')
    TT=Fnew;
    clear Fnew;
   case('salin');
    SS=Fnew;
   case('u-vel.');
    UU=Fnew;
   case('v-vel.');
    VV=Fnew;
  end
  
%keyboard
  f_chck=0;
  if f_chck>0
    [Fold,n,m,l] = read_hycom(finaOLD,finbOLD,vname);% read in 32-layer field
    Fold(Fold>0.0001*hg)=nan;
    IJ=[200,15;1020,15];
    sub_plot32v41(Fnew,Fold,HH,IJ,ZZn,ZZo,vname);
  end

  fprintf('  ####  Interpolate 1 fld: %9.2f sec\n\n',toc);

end

% 
fprintf('Saving interpolated files in mat files\n');
fmat=sprintf('%sTT_interp32to41.mat',pthout);
save(fmat,'TT');
fmat=sprintf('%sSS_interp32to41.mat',pthout);
save(fmat,'SS');
fmat=sprintf('%sUU_interp32to41.mat',pthout);
save(fmat,'UU');
fmat=sprintf('%sVV_interp32to41.mat',pthout);
save(fmat,'VV');

else  % f_mat
  fprintf('Loading interpolated files in mat files\n');
  fmat=sprintf('%sTT_interp32to41.mat',pthout);
  fprintf('%s\n',fmat);
  load(fmat);
  fmat=sprintf('%sSS_interp32to41.mat',pthout);
  fprintf('%s\n',fmat);
  load(fmat);
  fmat=sprintf('%sUU_interp32to41.mat',pthout);
  fprintf('%s\n',fmat);
  load(fmat);
  fmat=sprintf('%sVV_interp32to41.mat',pthout);
  fprintf('%s\n',fmat);
  load(fmat);
end;  

% ====================================
% Interpolation done
% Write output
%need to write 41 layers
%not 32
%flip through all variables and write by layers
%u-vel.
%v-vel.
%thknss
%temp
%salin
% ====================================
fprintf('Writing output ... \n');
% Output files:
hr2=3;
finaNEW = sprintf('%sarchv41lr_arcT%2.2i.%4.4i_%3.3i_%2.2i.a',...
		 pthout,TV,yr,iday,hr2);
finbNEW = sprintf('%sarchv41lr_arcT%2.2i.%4.4i_%3.3i_%2.2i.b',...
		 pthout,TV,yr,iday,hr2);

% Open files "relax" files (anom+mean old nest) for reading
% and nest files for writing
disp('Openning files ...');
fida_OLD = fopen(finaOLD,'r');  % 
fidb_OLD = fopen(finbOLD,'rt'); % use this to copy over *b
fida_NEW = fopen(finaNEW,'w');
fidb_NEW = fopen(finbNEW,'wt');


nlevU=0;  % level counter u-vel
nlevV=0;  % level counter v-vel
nlevD=0;  % level counter layer thickness 
nlevT=0;  % level counter layer T 
nlevS=0;  % level counter layer S 

if ~exist('TDENS','var')
%  TDENS=read_targ_dens(finbOLD,'salin','=',4);
  pthb='/home/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers/';
  flblk=sprintf('%sblkdat.input_ARCc0.08_41lev',pthb);
  TDENS=read_targ_dens_blkdat(flblk);
end

% Write header:
aa1='Daily nest from GLBb0.08 32-layer T07 interpolated into ARCc0.08 41 layers T11';


for nl=1:10
  aa=fgetl(fidb_OLD);
  if nl==4
    aa=aa1;
  end
  fprintf(fidb_NEW,[aa,'\n']);

  if nl==8
    dmm=aa(2:8);
    IDM=str2num(dmm);
  end
  if nl==9
    dmm=aa(2:8);
    JDM=str2num(dmm);
  end

  if exist('JDM','var')
  if (IDM~=nn | JDM~=mm)
    fprintf('ERR: Check dimensions in new grid  and old grid:\n');
    fprintf('New: IDM=%i, JDM=%i; Old: IDM=%i, JDM=%i\n',nn,mm,IDM,JDM);
    error('*** STOPPING ***');
  end
  end

end

IJDM=IDM*JDM;
npad=4096-mod(IJDM,4096);
toto=ones(npad,1);

% Write barotropic part that does not need changes:
while ischar(aa)
  aa=fgetl(fidb_OLD);

  if ~ischar(aa), break; end;

  I=strfind(aa,'=');
  bb=aa(I+1:end);
  S=sscanf(bb,'%f');
  bb=aa(1:I-1);
  wrd = deblank(sscanf(bb,'%s'));
  Imx = strmatch('v-vel.',aa);

  
% Copy over surface & barotropic fields
  if S(3)==0 | strmatch(wrd,'montg1')  %# hybr. layer #
    fprintf(' <= %s\n',aa);
% Read 1 record from HYCOM nest archv:    
    dmm  = fread(fida_OLD,IJDM,'float32','ieee-be');  % read field
    toto = fread(fida_OLD,npad,'float32','ieee-be');  % read npad 
% Write in the new nest file:
    fwrite(fida_NEW,dmm,'float32','ieee-be');
    fwrite(fida_NEW,toto,'float32','ieee-be');

    fprintf('==> %s\n',aa);
    fprintf(fidb_NEW,[aa,'\n']);
  end  % if S  - surface and barotropic fields are done
% either one:
  if strcmp(wrd,'u-vel.'); break; end;
  if ~isempty(Imx); break; end;
end

aa0=aa;

% Interpolated baroclinic part:
FLDS={'u-vel.';'v-vel.';'thknss';'temp';'salin'};
for kk=1:nlN
for ip=1:5
  wrd=FLDS{ip};
    
  Is=strfind(aa,'=');
  lns=length(wrd);
  aa=aa0;
  aa(1:Is-1)=' ';
  aa(1:lns)=wrd;
  
  fprintf('Field: %s, layer = %i\n',wrd,kk);
  
  switch (wrd)
% -----------------
%  Fields
% -----------------
   case('u-vel.')
    nlevU=nlevU+1;
    nlev=nlevU;
    sub_wfld2d_ab(nlev,UU,IDM,JDM,fida_NEW,fidb_NEW,aa,TDENS);
    
   case('v-vel.')
    nlevV=nlevV+1;
    nlev=nlevV;
    sub_wfld2d_ab(nlev,VV,IDM,JDM,fida_NEW,fidb_NEW,aa,TDENS);
    
   case('thknss')
    nlevD=nlevD+1;
    nlev=nlevD; 
    dmm=abs(squeeze(ZZn(nlev+1,:,:)-ZZn(nlev,:,:)));
    I=find(dmm<1e-3);
    dP=abs(dmm)*rg;
    dP(I)=nan;
    sub_wfld2d_ab(nlev,dP,IDM,JDM,fida_NEW,fidb_NEW,aa,TDENS);
    
   case('salin')
    nlevS=nlevS+1;
    nlev=nlevS;
    sub_wfld2d_ab(nlev,SS,IDM,JDM,fida_NEW,fidb_NEW,aa,TDENS);

   case('temp')
    nlevT=nlevT+1;
    nlev=nlevT;
    sub_wfld2d_ab(nlev,TT,IDM,JDM,fida_NEW,fidb_NEW,aa,TDENS);
  end % switch wrd - field name 	 

end;  % for ip variables  
end % reading *b file, ischar aa

fclose(fida_NEW);
fclose(fidb_NEW);
fclose(fida_OLD);
fclose(fidb_OLD);

fprintf(' OUTPUT saved: %s\n',finaNEW);
fprintf(' OUTPUT saved: %s\n\n',finbNEW);





