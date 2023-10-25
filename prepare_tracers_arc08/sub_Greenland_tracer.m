function Rmask = sub_Greenland_tracer(BB);
% For multiple tracers
% write each tracer by month
% for ktr=1:ntracr
% ntracr is specified in blkdat.input
ktrcr    = BB.Tracer;  % current tracer #
PTH      = BB.PTH;
ftrca    = BB.output_afile;
ftrcb    = BB.output_bfile;
trca_fid = BB.output_afid;
trcb_fid = BB.output_bfid;
HH       = BB.hycom_topo;
elon     = BB.hycom_LON;
alat     = BB.hyom_LAT;
cf       = BB.area_cf;       % 16e6 - mean cell area, m2, for runoff m/s->m3/sec
ACell    = BB.Area_cell_m2;  % if empty, - use area_cf;		      
TDENS    = BB.TDENS;
Rmask    = BB.Rmask;
ntopo    = BB.TopoN;
iyr      = BB.Year_River;

if ~isempty(ACell),
  cf=[];
end

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
nlev=length(TDENS);
%nlev=0;
%if nlev==0
%  fss=sprintf('%srelax_sal.b',PTH.data);
%  TDENS=read_targ_dens(fss,'sal','=',3);
%  nlev=length(TDENS);
%end


strB=' trc: month,layer,dens,range = ';



% Rivers:
flriva=sprintf('%srivers_%2.2i_Greenland_%i.a',PTH.river,ntopo,iyr);
flrivb=sprintf('%srivers_%2.2i_Greenland_%i.b',PTH.river,ntopo,iyr);
riv_fid=fopen(flriva,'r','ieee-be');
if ~exist(flriva,'file')
  fprintf('River files not found: %s\n',flriva);
  error('Check file');
end


% Find rivers that freeze up during year
% Use 10% of annual mean values during
% freezing time
% Problem: 0 tracer will wipe out other tracers
% if they end up in the frozen domain
% forced to relax to 0
[Rmean,JR] = sub_rivermean(riv_fid,GRMSK);
IALL=find(Rmean>0);
Rmask(IALL) = 1;

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
    error('# of river points is in imbalance ...');
  end
  
  if ~isempty(IB)
    A(IB)=0.1*Rmean(IB);
  end
  J=find(~isnan(A));
  
%keyboard
%  CH=[];
%  kch=0;
  if ~isempty(cf)
    flx=A(J)*cf;
  else
    flx=A(J).*ACell(J); %m3/sec - "flux" of tracer
  end
  
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

  chck_plt=0;
  if chck_plt>0
% Subsample Greenland
    T=TRC(jg1:jg2,ig1:ig2);
    H=HH(jg1:jg2,ig1:ig2);
%    np=NP(jg1:jg2,ig1:ig2);
    F=A(jg1:jg2,ig1:ig2);
%    Ts=As(jg1:jg2,ig1:ig2);
    T(T==0)=nan;
    clf;
    contour(H,[0 0],'Color',[0.9 0.9 0.9]);
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
  fprintf('Tracer %i:  Writing tracers %s\n',ktrcr,ftrcb);
  TRC(isnan(TRC))=0;
  TRC=TRC';
  Btot=reshape(TRC,IJDM,1);
  fprintf('Tracer %i:  Writing tracers %s\n',ktrcr,ftrca);
  nmix=2;  % over what depth distribute tracers
  for kl=1:nlev
    if kl>nmix, 
      B=Btot*0; 
    else
      B=Btot/nmix;
    end;
    
    btxt=sprintf('%s%2.2i  %2.2i %6.3f %14.7e %14.7e\n',...
		 strB,k,kl,TDENS(kl),min(B),max(B));
    fprintf(trcb_fid,btxt);
    
    fwrite(trca_fid,B,'float32');
    fwrite(trca_fid,toto,'float32'); % padding at the end of record
  end;
  
  
end;
  
fclose(riv_fid);



return