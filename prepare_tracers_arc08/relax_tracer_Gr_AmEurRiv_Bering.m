% Tracers in hycom: 4 or 5 tracers - track 
% 1)Greenland runoff and
% 2) AMerican/
% 3) Eurasian river particles: Can be splitted into:
%    (3) Eastern Euras. Rivers (Kolyma, Khatanga, Lena)
%    (4) Western Euras. Rivers (Yenisei, Ob, Pechora, S.Dvina)
% 4) (or (5) if Euras. Riv. splitted)  Bering strait 
% seeded at river points specified in 
% HYCOM river files
%
% River fluxes from HYCOM are converted from m/s -> m3/s
% These values are set as max conc. value for tracers
% for Bering Str.: max conc = 1 for all points across the strait
%
%
% Multiple monthly tracer climatology 
% should be written in the following order:
% for month=1,12
%  for tracers = 1,ntracr
%   for layers = 1,nlev
%
%  see forfun.f: subroutine rdrlax
%  Skipping not needed months (skip forward to desired month)
%  do ktr=1,ntracr
%    do k=1,kk  <-- # of vert. layers
%      call skmonth(iunit) <--- skips monthly 2D field
%
% Concentration made proportional to river
% flux at the grid
% distribute particles within the layer
% 'thkriv' - nominal thickness of river inflow specified
% in blkdat.input = 6m in this experiment
% ~ 2 upper layers - Note in 41-layer version
%       upper layers are thinner (~ 1m thick)
% below - all tracers are 0
%
% Approach: find river grid cells
% along North-East Gr. coast
% using river input *.ab files
% (prepared from Bamber data)
% Put tracers in the same grid points
% the concentration is proportional 
% to the river runoff
%
% NOTE: this is not the tracer flux in the model
% HYCOM relaxes tracer concentration to these values!
% if tr. value in the ocean grid cell = relaxes value,
%  tracer flux = 0, no relaxation needs to be done
% 
% Tracer relaxation file is actual "concentration", unitless
% this is NOT A FLUX 
% to get the actual tracer flux:
% at every tracer input point for dayN (C2) with tracer relax
% and same day without tracer relax (C1)
% tracer conc. change: dC=C2-C1
% Flux = dC/dTime*Volume_grid_cell [kg/sec]
%
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close


%s = matlabpool('size');
%if s==0
%  matlabpool open 12;
%end


regn='ARCc0.08';
expt='100';  % hycom experiment *10
ntopo=11;      % HYCOM topo version
YearGr=2005;  % Gr. river runoff river


cf=16e6; % to get roughly m3/sec integrated over 1grid cell
Nrlx=25;  % # of relax. grid cells

hg=2^100; 

PTH.data='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/force/relax/';
PTH.topo='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
PTH.river='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/force/rivers/';

mday=[31;28;31;30;31;30;31;31;30;31;30;31];
if mod(YearGr,4)==0,
  mday(2)=29;
end

% Get topo and grid:
fltopo=sprintf('%sdepth_ARCc0.08_%2.2i.nc',PTH.topo,ntopo);
HH   = nc_varget(fltopo,'Bathymetry');
alat = nc_varget(fltopo,'Latitude');
elon = nc_varget(fltopo,'Longitude');
LAT  = alat;
LON  = elon;
[m,n]= size(HH);
[mm,nn]= size(HH);
IDM=n;
JDM=m;
IJDM=IDM*JDM;
npad=4096-mod(IJDM,4096);
toto=ones(npad,1);


%DX=zeros(mm,nn);
%DY=zeros(mm,nn);
%for i=1:nn-1
%  dx=distance_spheric_coord(LAT(:,i),LON(:,i),LAT(:,i+1),LON(:,i+1));
%  DX(:,i)=dx;
%end
%DX(:,nn)=dx;
%for j=1:mm-1
%  dy=distance_spheric_coord(LAT(j,:),LON(j,:),LAT(j+1,:),LON(j+1,:));
%  DY(j,:)=dy;
%end
%DY(mm,:)=dy;

[DX,DY] = sub_dx_dy(LON,LAT);
ACell=DX.*DY;
%keyboard

% Read t. densities 41 layers:
%fss=sprintf('%srelax_sal.b',PTH.data); % 32 layers
%TDENS=read_targ_dens(fss,'sal','=',3);
if ~exist('TDENS','var')
%  TDENS=read_targ_dens(finbOLD,'salin','=',4);
  pthb='/home/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers/';
  flblk=sprintf('%sblkdat.input_ARCc0.08_41lev',pthb);
  TDENS=read_targ_dens_blkdat(flblk);
end

nlev=length(TDENS);
fprintf('# v. levels in HYCOM=%i\n',nlev);

% ================
% Tracer files:
% OPen files to write tracer conc.
ftrca = sprintf('%srelax_trcr_dummy.a',PTH.data);
ftrcb = sprintf('%srelax_trcr_dummy.b',PTH.data);
trca_fid=fopen(ftrca,'w','ieee-be');
trcb_fid=fopen(ftrcb,'wt');
ftrcaDUM = ftrca;
ftrcbDUM = ftrcb;

% Write heading in *.b:
fprintf('Writing %s\n',ftrcb);
%aa=sprintf(...
%  'Passive tracers Greenland, Amer.Riv., Euras.Riv., Bering');
aa=sprintf(...
  'Passive tracers Greenland, Amer.Riv., East/West Euras.Riv., Bering');
disp(aa);
fprintf(trcb_fid,[aa,'\n']);

aa = sprintf('Topo=%2.2i, N.Vert Lev=%i',ntopo,nlev);
disp(aa);
fprintf(trcb_fid,[aa,'\n']);

aa='Tracers seeded in river grid cells';
disp(aa);
fprintf(trcb_fid,[aa,'\n']);

aa='Trac.Conc = runoff m3/sec intgr over 1 grid cell';
disp(aa);
fprintf(trcb_fid,[aa,'\n']);

%aa='# of grid cells where tracers seeded is same as hycom river';
%disp(aa);
%fprintf(trcb_fid,[aa,'\n']);

aa=sprintf('i/jdm = %i %i',IDM,JDM);
disp(aa);
fprintf(trcb_fid,[aa,'\n']);
% ================

% First, in a temporary file
% Tracers are written into relax.trcr.a, b
% by tracers: layers 1:kk, month=1:12
%  Then file is rearranged in the right order (Months, layers, tracers)
%keyboard
Rmask = zeros(mm,nn);
% ============================
% Tracer # 1 - Greenland Tracer
% ============================
BB.Tracer       = 1;
BB.PTH          = PTH;
BB.output_afile = ftrca;
BB.output_bfile = ftrcb;
BB.output_afid  = trca_fid;
BB.output_bfid  = trcb_fid;
BB.TopoN        = ntopo;
BB.hycom_topo   = HH;
BB.hycom_LON    = elon;
BB.hyom_LAT     = alat;
BB.area_cf      = 16e6;  % 16e6 - mean cell area, m2, for runoff m/s->m3/sec
BB.Area_cell_m2 = ACell; % if empty, - use area_cf;
BB.TDENS        = TDENS;
BB.Year_River   = YearGr; % Yr for Gr. rivers
BB.Rmask        = Rmask; % river sources mask

Rmask = sub_Greenland_tracer(BB);

% ============================
% Tracer # 2 - American River (Mackenzie):
% ============================
BB.Tracer       = 2;
BB.river_flag   = 'amer';
BB.Rmask        = Rmask; % river sources mask
Rmask = sub_AmerEurRiv_tracer(BB);
%keyboard
% ============================
% Tracer # 3 - Eurasian Rivers (S.Dvina, Pechora
%              Yenisey, Pyasina, Olenek, Lena,
%              Kolyma):
% ============================
%BB.Tracer       = 3;
%BB.river_flag   = 'euras';
%BB.Rmask        = Rmask; % river sources mask
%Rmask = sub_AmerEurRiv_tracer(BB);

% ============================
% Tracer # 3 - East Eurasian Rivers (Lena,Khatanga,
%              Kolyma):
% ============================
BB.Tracer       = 3;
BB.river_flag   = 'E-euras';
BB.Rmask        = Rmask; % river sources mask
Rmask = sub_AmerEurRiv_tracer(BB);

% ============================
% Tracer # 4 - West Eurasian Rivers (Kara/Bar.Seas:
%              Yenisey, Ob, Sev. Dvina, Pechora)
% ============================
BB.Tracer       = 4;
BB.river_flag   = 'W-euras';
BB.Rmask        = Rmask; % river sources mask
Rmask = sub_AmerEurRiv_tracer(BB);

% ============================
% Tracer # 5 - Bering Strait:
% ============================
BB.Tracer       = 5;
BB.river_flag   = 'bering';
BB.Rmask        = Rmask; % river sources mask
Rmask = sub_Bering_tracer(BB);

ntracr=BB.Tracer;
fprintf('# of tracers: %i\n',ntracr);

% Close tracer relax files
fclose(trca_fid);
fclose(trcb_fid);

fprintf(' \nDummy Tracer files are ready \n\n');


% ========================
% Rearrange HYCOM files:
% Dummy: for itr=1:ntracr (for im=1:12 (for k=1:nlev,...))
% to:    for im=1:12 (for itr=1:ntracr (for k=1:nlev ...)) 
ftrcaOLD = ftrcaDUM;
ftrcbOLD = ftrcbDUM;
fidaOLD  = fopen(ftrcaOLD,'r','ieee-be');
fidbOLD  = fopen(ftrcbOLD,'rt');

ftrca    = sprintf('%srelax_trcr%2.2i_Greenl_Riv_Bering_T%2.2i.a',...
		   PTH.data,ntracr,ntopo);
ftrcb    = sprintf('%srelax_trcr%2.2i_Greenl_Riv_Bering_T%2.2i.b',...
		   PTH.data,ntracr,ntopo);
fidaNEW  = fopen(ftrca,'w','ieee-be');
fidbNEW  = fopen(ftrcb,'wt');

% Copy header:
for nl=1:5
  aa=fgetl(fidbOLD);
  fprintf(fidbNEW,[aa,'\n']);
end

strB=' trc: month,layer,dens,range = ';

%keyboard


for im=1:12
  fprintf('\nMonth %i\n',im);
  for ktr=1:ntracr
    for k=1:nlev
      fprintf('Reading/writing mo=%i, tracer=%i, lev=%i\n',im,ktr,k);
      nrec=(ktr-1)*(12*nlev)+(im-1)*nlev+k-1; % this should match field order in dummy
%      nrec=(im-1)*nlev*ntracr+(ktr-1)*nlev+k-1; % this should match field order in dummy
      stat=fseek(fidaOLD,nrec*(IJDM+npad)*4,-1);
      if stat<0,
	[msg,errnm]=ferror(fidaOLD);
	fprintf('%s\n',msg);
	error('Couldnot find record location in %s\n',ftrcaOLD);
      end
      
      B=fread(fidaOLD,IJDM,'float32');  % read 2D field (1 layer)
%      dm1=fread(fidaOLD,npad,'float32','ieee-be');  % read npad 
      fwrite(fidaNEW,B,'float32');
      fwrite(fidaNEW,toto,'float32'); % padding at the end of record

      btxt=sprintf('%s%2.2i  %2.2i %6.3f %14.7e %14.7e\n',...
		 strB,im,k,TDENS(k),min(B),max(B));
      fprintf(fidbNEW,btxt);
    end
  end
end

% Close tracer relax files
fclose(fidaNEW);
fclose(fidbNEW);
fclose(fidaOLD);
fclose(fidbOLD);

fprintf(' \nTracer files are ready %s \n\n',ftrca);

%keyboard

scmd=sprintf('rm %s',ftrcaOLD);
system(scmd);
scmd=sprintf('rm %s',ftrcbOLD);
system(scmd);


% 
% ===================================
% Prepare relaxation time files:
fina = sprintf('%srelax_rmutr_Greenl_Riv_Bering_T%2.2i.a',PTH.data,ntopo);
finb = sprintf('%srelax_rmutr_Greenl_Riv_Bering_T%2.2i.b',PTH.data,ntopo);
fina_fid = fopen(fina,'w','ieee-be');
finb_fid = fopen(finb,'wt');


% Relaxation mask
% here - 1 mask for all tracers
% can be several masks, 
% then need to change 1 line in relax_rmutr.b
% to "Relaxation Masks" (from "Relaxation Mask")
% Relaxation time scale
% is 1/1day in the river grid cells
rmu1=1./86400; %  s^-1
%rmu2=1./(5*86400); %


RMU=Rmask*0;
RMU(Rmask>0)=rmu1;
%RMU(JR)=rmu2;
%r=RMU(jg1:jg2,ig1:ig2);
Ir=find(RMU>0);
nIr=length(Ir);

fprintf(finb_fid,'Relaxation Mask, Tracers Greenl, Amer. Euras.Rivs., Bering\n');
fprintf(finb_fid,'Total # river points, %i points \n',nIr);
fprintf(finb_fid,' \n\n');
fprintf(finb_fid,'i/jdm = %i  %i\n',IDM,JDM);
fprintf(finb_fid,' rmutr: range =  %15.7e %15.7e\n',...
      min(min(RMU)), max(max(RMU)) );

fclose(finb_fid);

%r=RMU(jg1:jg2,ig1:ig2);
%pcolor(r); shading flat;


RMU=RMU';
A=reshape(RMU,IJDM,1);
fwrite(fina_fid,A,'float32','ieee-be');
fwrite(fina_fid,toto,'float32','ieee-be');  % padding at the end
fclose(fina_fid);

fprintf(' \nTracer relaxation time files are ready \n');
fprintf('%s\n',fina);

%fsv=0;
%if fsv>0
%  fmat=sprintf('%srelax_tracer.mat',PTH.data);
%  fprintf('saving %s\n',fmat);
%  save(fmat,'NP');
%elseif fsv==2
%  load(fmat);
%end

%matlabpool close;








