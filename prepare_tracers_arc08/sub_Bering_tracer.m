function Rmask = sub_Bering_tracer(BB);
% Bering tracers:
% The total FW flux through Bering strait
% is ~2500 km3/yr (0.08Sv) - Woodgate & Aagaard, 2005
% up to ~3200 km3/yr (0.1Sv) Woodgate et al., 2015?
%
% It is perhaps better to use the total volume flux 
% and redistribute it over the grid cells where tracers are intorduces
%
% Here, tracer concentration is set = 1
% in all grid points in all layers

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
f_riv    = BB.river_flag;  % specify rivers: 'amer' or 'euras'
Rmask    = BB.Rmask;

if ~isempty(ACell),
  cf=[];
end

[m,n]= size(HH);
[mm,nn]= size(HH);
nlev=length(TDENS);

% Bering strait indices:
% a box where tracer is released
i1=639;
i2=656;
j1=1919;
j2=1925;

IDM=n;
JDM=m;
IJDM=IDM*JDM;
npad=4096-mod(IJDM,4096);
toto=ones(npad,1);


strB=' trc: month,layer,dens,range = ';

fprintf('Tracers Bering Strait \n');

TRC=zeros(mm,nn);
TRC(j1:j2,i1:i2)=1.0;
TRC(HH>=0)=0;
Irv=find(TRC>0);
Rmask(Irv)=1;
TRC=TRC';
B=reshape(TRC,IJDM,1);

for k=1:12  % time - 12mo
  fprintf('sub_Bering_tracer:    Reading month %i\n',k);

  fprintf('Bering: Writing tracers %s\n',ftrcb);
  fprintf('Bering: Writing tracers %s\n',ftrca);

  for kl=1:nlev
    btxt=sprintf('%s%2.2i  %2.2i %6.3f %14.7e %14.7e\n',...
		 strB,k,kl,TDENS(kl),min(B),max(B));
    fprintf(trcb_fid,btxt);
    
    fwrite(trca_fid,B,'float32');
    fwrite(trca_fid,toto,'float32'); % padding at the end of record
  end;
  
  
end;
  

return