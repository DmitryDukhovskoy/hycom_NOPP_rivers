function Rmask = sub_AmerEurRiv_tracer(BB);
% Write American River (Mackenzie) Tracer
% or Eurasian Rivers (Dvina,Pechora,Ob,Yenisei <- missing
%  Pyasina, Olenek, Lena, Kolyma) Tracers
%
ktrcr    = BB.Tracer;  % current tracer #
PTH      = BB.PTH;
ftrca    = BB.output_afile;
ftrcb    = BB.output_bfile;
trca_fid = BB.output_afid;
trcb_fid = BB.output_bfid;
HH       = BB.hycom_topo;
elon     = BB.hycom_LON;
alat     = BB.hyom_LAT;
ntopo    = BB.TopoN;
cf       = BB.area_cf;       % 16e6 - mean cell area, m2, for runoff m/s->m3/sec
ACell    = BB.Area_cell_m2;  % if empty, - use area_cf;		      
TDENS    = BB.TDENS;
f_riv    = BB.river_flag;  % specify rivers: 'amer' or 'euras'
Rmask    = BB.Rmask;

if ~isempty(ACell),
  cf=[];
end

[m,n]= size(HH);
nlev=length(TDENS);

% Reading rivers & bathymetry:
IDM=n;
JDM=m;
IJDM=IDM*JDM;
npad=4096-mod(IJDM,4096);
toto=ones(npad,1);


%flriva=sprintf('%srivers_09_corrected.a',PTH.river);
%flrivb=sprintf('%srivers_09_corrected.b',PTH.river);
flriva=sprintf('%srivers_%2.2i.a',PTH.river,ntopo);
flrivb=sprintf('%srivers_%2.2i.b',PTH.river,ntopo);
%riv_fid=fopen(flriva,'r');
riv_fid=fopen(flriva,'r','ieee-be');

strB=' trc: month,layer,dens,range = ';

RMSK=zeros(JDM,IDM);

if strncmp(f_riv,'amer',4)
  RMSK(1600:1700,400:450)=1;  % Mackenzie river
elseif strncmp(f_riv,'euras',4) % All Eur. rivers 1 tracer
  RMSK(1752:1904,781:1108)=1; % E-Sib. sea
  RMSK(1610:1795,1000:1315)=1; % E. Lapt. sea
  RMSK(640:1609,1226:1597)=1; % all other seas
  RMSK(640:1322,1220:1355)=0; % notneeded area
elseif strncmp(f_riv,'E-euras',5) % East Euras. rivers
  RMSK(1752:1904,781:1108)=1; % E-Sib. sea
  RMSK(1610:1795,1000:1315)=1; % E. Lapt. sea
  RMSK(1452:1609,1226:1406)=1; % W. Lapt. sea
elseif strncmp(f_riv,'W-euras',5) % West Euras. rivers
  RMSK(640:1451,1226:1598)=1; % Barents, Kara Seas
  RMSK(810:1322,1220:1355)=0; % notneeded area
else
  error('sub_AmerEurRiv_tracer.m:  Unknown river flag  ...');
end

fprintf('Tracers for %s rivers \n',f_riv);

for k=1:12  % time - 12mo
  fprintf('sub_AmerEurRiv_tracer:    Reading month %i\n',k);

  A  = fread(riv_fid,IJDM,'float32'); % read 2D field
  dm1= fread(riv_fid,npad,'float32');  % Padding = size(toto)
  if size(dm1) ~= size(toto)
    error('Padding in HYCOM file ???');
  end
  
  I=find(A>1e10);
  A(I)=NaN;
  A=reshape(A,IDM,JDM)';
  A(A==0)=nan;
  A(RMSK==0)=nan;

  TRC=HH*0;
  J=find(~isnan(A));
  
%  flx=A(J)*cf;
%  TRC(J)=round(flx);
  if ~isempty(cf)
    flx=A(J)*cf;
  else
    flx=A(J).*ACell(J); % m3/sec
  end
  
  TRC(J)=flx; % trac.flux per 1 gr.cell over whole depth
  Irv=find(TRC>0); % mask river points
  Rmask(Irv)=1;
    
% Check that none of river/tracer sources are onland
  TRC(TRC==0)=nan;
  JT=find(~isnan(TRC));
  IL=find(HH(JT)>=0);
  if ~isempty(IL),
    error('sub_AmerEurRiv_tracer.m: Need to move river source off land');
  end
  
  fprintf('AmerEurRiv: Writing tracers %s\n',ftrcb);
  TRC(isnan(TRC))=0;
  TRC=TRC';
  Btot=reshape(TRC,IJDM,1);
  fprintf('AmerEurRiv: Writing tracers %s\n',ftrca);
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
  

% Rivers:
%riv_fid=fopen(flriva,'r');
%fseek(riv_fid,6*4*(npad+IJDM),-1);
%dmm=fread(riv_fid,1,'int','ieee-be');
%[A,counta]=fread(riv_fid,IJDM,'float32','ieee-be');
%I=find(A>1e10);
%A(I)=NaN;
%A=reshape(A,IDM,JDM)';



%keyboard

chck=0;
if chck>0
  dmm=Btot;
  dmm=reshape(dmm,[IDM,JDM])';
  A=dmm./nmix; % tracer flux. in 1 gr. cell, per layer, m3/s
  tmp=log10(A);
  tmp(A==0)=NaN;
  
  figure(10);
  clf
  contour(HH,[0 0],'g');
  hold on;
  pcolor(tmp); shading flat;
  caxis([-10 -2]);
  colorbar
  title('log10(River runoff), m/s')
  keyboard
end



return