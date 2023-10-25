function RR = river_FWflux_byTrcrs(BB);
% Create River runoff and Bering FW fluxes
% m3/mo for plotting tracer concentration
% in terms of FW content and S change
% FW fluxes from monthly climatology
PTH.data='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/force/relax/';
PTH.topo='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
PTH.river='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/force/rivers/';

HH       = BB.hycom_topo;
elon     = BB.hycom_LON;
alat     = BB.hyom_LAT;
f_riv    = BB.river_flag;  % specify rivers: 'amer' or 'euras'
ACell    = BB.Area_cell_m2;  % if empty, - use area_cf;	
yr1      = 1993;
yr2      = 2016;
ntopo    = 11;

mday = [31;28;31;30;31;30;31;31;30;31;30;31];

[m,n]= size(HH);
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

fprintf('FWFlux %s rivers \n',f_riv);

clear FWF
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

  J=find(~isnan(A));
  
%  flx=A(J)*cf;
%  TRC(J)=round(flx);
  flx=A(J).*ACell(J); % m3/sec
  FWF(k)=nansum(flx)*3600*24*mday(k); % m3/mo
  
%  if k==7, keyboard; end;
end;
  
fclose(riv_fid);

cc=0;
for iyr=yr1:yr2
  for im=1:12
    cc=cc+1;
    RR.FWF_m3_mo(cc,1)=FWF(im);
    RR.TM(cc,1)=datenum(iyr,im,15);
  end
end


return