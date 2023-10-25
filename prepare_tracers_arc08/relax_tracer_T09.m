% Tracers in hycom - track Greenland river particles
% seeded at river points specified in 
% HYCOM river files
% Concentration make proportional to river
% flux at the grid
% distribute particles within the layer
% 'thkriv' - nominal thickness of river inflow specified
% in blkdat.input = 6m in this experiment
% ~ 2 upper layers
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
expt='050';  % hycom experiment *10
ntopo=9;      % HYCOM topo version
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
fltopo=sprintf('%sdepth_ARCc0.08_09.nc',PTH.topo);
HH   = nc_varget(fltopo,'Bathymetry');
alat = nc_varget(fltopo,'Latitude');
elon = nc_varget(fltopo,'Longitude');
[m,n]= size(HH);

% Greenland subdomain:
jg1=320;
jg2=1120;
ig1=500;
ig2=1000;


% HYCOM dim:
IDM=n;
JDM=m;
IJDM=IDM*JDM;
npad=4096-mod(IJDM,4096);
toto=ones(npad,1);


% Put tracers along Greenland 
DL1=find(elon>-5);
DL2=find(elon<-75);
DL3=find(alat<55);
DL4=find(alat>85);
GRMSK=ones(JDM,IDM);
GRMSK(DL1)=0;
GRMSK(DL2)=0;
GRMSK(DL3)=0;
GRMSK(DL4)=0;
GRMSK(:,1:500)=0;
GRMSK(:,980:end)=0;


% # of v. levels - read it 
% specify
nlev=0;
if nlev==0
  fss=sprintf('%srelax_sal.b',PTH.data);
  TDENS=read_targ_dens(fss,'sal','=',3);
  nlev=length(TDENS);
end

% OPen files to write tracer conc.
ftrca = sprintf('%srelax_trcr_Greenland.a',PTH.data);
ftrcb = sprintf('%srelax_trcr_Greenland.b',PTH.data);
trca_fid=fopen(ftrca,'w','ieee-be');
trcb_fid=fopen(ftrcb,'wt');

% Write heading in *.b:
fprintf('Writing %s\n',ftrcb);
aa=sprintf(...
  'Passive tracers Greenl coast, Bamber Gr. runoff %i',YearGr);
disp(aa);
fprintf(trcb_fid,[aa,'\n']);

aa='Tracers seeded in river grid cells';
disp(aa);
fprintf(trcb_fid,[aa,'\n']);

aa='Concentration is runoff m3/sec intgr over 1 grid cell';
disp(aa);
fprintf(trcb_fid,[aa,'\n']);

aa='# of grid cells where tracers seeded is same as hycom river';
disp(aa);
fprintf(trcb_fid,[aa,'\n']);

aa=sprintf('i/jdm = %i %i',IDM,JDM);
disp(aa);
fprintf(trcb_fid,[aa,'\n']);

strB=' trc: month,layer,dens,range = ';



% Rivers:
flriva=sprintf('%srivers_09_Greenland_2004.a',PTH.river);
flrivb=sprintf('%srivers_09_Greenland_2004.b',PTH.river);
riv_fid=fopen(flriva,'r','ieee-be');

% Find rivers that freeze up during year
% Use 10% of annual mean values during
% freezing time
% Problem: 0 tracer will wipe out other tracers
% if they end up in the frozen domain
% forced to relax to 0
[Rmean,JR] = sub_rivermean(riv_fid,GRMSK);
IALL=find(Rmean>0);

chck=0;
if chck>0
  clf;
  contour(HH,[0 0],'c');
  hold
  dmm=Rmean;
  dmm(JR)=nan;
  In=find(dmm>0);
  for k=1:length(JR);
    [j,i]=ind2sub([JDM,IDM],JR(k));
    plot(i,j,'m*');
  end
  for k=1:length(In);
    [j,i]=ind2sub([JDM,IDM],In(k));
    plot(i,j,'ko');
  end
  set(gca,'xlim',[ig1 ig2],'ylim',[jg1 jg2]);
end;

frewind(riv_fid);
%
% Some rivers freeze during winter
%NP=HH*0;  % keep track of grid cells with tracers
for k=1:12  % time - 12mo
  fprintf('Reading month %i\n',k);

  A  = fread(riv_fid,IJDM,'float32'); % read 2D field
  dm1= fread(riv_fid,npad,'float32');  % Padding = size(toto)
  if size(dm1) ~= size(toto)
    error('Padding in HYCOM file ???');
  end
  
  I=find(A>1e10);
  A(I)=NaN;
  A=reshape(A,IDM,JDM)';
  A(A==0)=nan;
  A(GRMSK==0)=nan;
  
%  if k<7; continue; end
%keyboard
  
  TRC=HH*0;
  J=find(~isnan(A));
  j0=899;
  i0=900;

% Fix frozen rivers:  
  nR=length(J)/length(IALL)*100;
  fprintf('   %7.3f%% rivers active\n',nR);
  dmm=Rmean;
  dmm(J)=nan;
  dmm(dmm==0)=nan;
  IB=find(~isnan(dmm));
  mtt=length(IB)+length(J);
  if mtt~=length(IALL),
    error('# of river points disbalance ...');
  end
  
  if ~isempty(IB)
    A(IB)=0.1*Rmean(IB);
  end
  J=find(~isnan(A));
  

%  CH=[];
%  kch=0;
  flx=A(J)*cf;
  TRC(J)=flx;
%  parfor ip=1:size(J,1);
%    if mod(ip,round(0.25*size(J,1)))==0
%%      fprintf('== River  mask %5.1f%%...\n',ip/size(J,1)*100);
%      fprintf('== Parallel loop: River  mask, ip %i...\n',ip);
%    end
%    
%    [j,i]=ind2sub(size(A),J(ip));
%    flx=A(j,i)*cf;  % ~m3/sec for the grid cell
%    DMM(ip)=DMM(ip)+flx;
%    TRC(j,i)=TRC(j,i)+flx;
%  end
%  TRC(J)=DMM;
  TRC(HH>=0)=nan;
%  NP(J)=NP(J)+1;

  chck_plt=1;
  if chck_plt>0
% Subsample Greenland
    T=TRC(jg1:jg2,ig1:ig2);
    H=HH(jg1:jg2,ig1:ig2);
%    np=NP(jg1:jg2,ig1:ig2);
    F=A(jg1:jg2,ig1:ig2);
%    Ts=As(jg1:jg2,ig1:ig2);
    T(T==0)=nan;
    clf;
    contour(H,[0 0],'k');
    hold on;
    pcolor(T); shading flat;
    set(gca,'xlim',[50 420],'ylim',[100 750]);
    caxis([0 200]);
    colorbar;
    spp=sprintf('Tracer conc. map, mo=%i', k);
    title(spp);
    drawnow
  end  
  
%  keyboard
  
  fprintf(' Writing tracers %s\n',ftrcb);
  TRC(isnan(TRC))=0;
  TRC=TRC';
  B=reshape(TRC,IJDM,1);
  fprintf(' Writing tracers %s\n',ftrca);
  for kl=1:nlev
    if kl>2, B=B*0; end;
    btxt=sprintf('%s%2.2i  %2.2i %6.3f %14.7e %14.7e\n',...
		 strB,k,kl,TDENS(kl),min(B),max(B));
    fprintf(trcb_fid,btxt);
    
    fwrite(trca_fid,B,'float32');
    fwrite(trca_fid,toto,'float32'); % padding at the end of record
  end;
  
  
end;
  
fclose(riv_fid);
fclose(trca_fid);
fclose(trcb_fid);

fprintf(' \nTracer files are ready \n\n');

% 
% Prepare relaxation time files:
fina = sprintf('%srelax_rmutr_Greenland.a',PTH.data);
finb = sprintf('%srelax_rmutr_Greenland.b',PTH.data);
fina_fid = fopen(fina,'w','ieee-be');
finb_fid = fopen(finb,'wt');

%
% Relaxation time scale
% is 1/1day in the river grid cells
% that are permenet
% and weaker over the grid cells
% where river is seasonal
rmu1=1./86400; %  s^-1
rmu2=1./(5*86400); %

% Mask of tracer points:
[a1,a2]=size(Rmean);
if a1~=JDM,
  error('Check size Rmean ...');
end

RMU=Rmean*0;
RMU(Rmean>0)=rmu1;
RMU(JR)=rmu2;
%r=RMU(jg1:jg2,ig1:ig2);
Ir=find(RMU>0);
nIr=length(Ir);

fprintf(finb_fid,'Relaxation Mask, Tracers Greenland\n');
fprintf(finb_fid,'Along Greenland Coast, %i points \n',nIr);
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

fprintf(' \nTracer relaxation time files are ready \n\n');

%fsv=0;
%if fsv>0
%  fmat=sprintf('%srelax_tracer.mat',PTH.data);
%  fprintf('saving %s\n',fmat);
%  save(fmat,'NP');
%elseif fsv==2
%  load(fmat);
%end

%matlabpool close;








