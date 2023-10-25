% Tracers in hycom - track Greenland river particles
% Check hycom binary prepared in relax_tracer.m
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close

regn='ARCc0.08';
ntopo=9;      % HYCOM topo version
kplt=1; % layer to plot



PTH.data='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/force/relax/';
PTH.topo='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
PTH.river='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/force/rivers/';

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

ftrca = sprintf('%srelax_trcr_Greenland.a',PTH.data);
ftrcb = sprintf('%srelax_trcr_Greenland.b',PTH.data);

fidb = fopen(ftrcb,'r');  % read I,J from *.b
for nl=1:5
  aa=fgetl(fidb);
  disp(aa);
end

is=strfind(aa,'= ');
[ID,JD] = strread(aa(is+1:end),'%d%d');

disp(['Grid I=',num2str(ID),' J=',num2str(JD)]);
IJDM=ID*JD;
npad=4096-mod(IJDM,4096);
toto=ones(npad,1);

% # of v. levels - read it 
% specify
nlev=0;
if nlev==0
  fss=sprintf('%srelax_sal.b',PTH.data);
  TDENS=read_targ_dens(ftrcb,'trc','=',3);
  nlev=length(TDENS);
end

fida = fopen(ftrca,'r');  % 
for im=12:12
  fprintf('Reading mo=%i lev=%i\n',im,kplt);
  nrec=(im-1)*nlev+kplt-1;
  stat=fseek(fida,nrec*(IJDM+npad)*4,-1);
  dmm=fread(fida,IJDM,'float32','ieee-be');  % read 2D field (1 layer)
%  dm1=fread(fid1,npad,'float32','ieee-be');  % read npad 

  A=reshape(dmm,ID,JD);
  A=A';
  
  ll=find(A>0);
  fprintf('# of river cells %i\n',length(ll));

  mnA=min(min(A(A>0)));
  mxA=max(max(A));
  
  H=HH(jg1:jg2,ig1:ig2);
  F=A(jg1:jg2,ig1:ig2);
%    Ts=As(jg1:jg2,ig1:ig2);
  F(F==0)=nan;
  
  dmm=A;
  dmm(jg1:jg2,ig1:ig2)=nan;
  I=find(dmm>0);
  if ~isempty(I),
    fprintf('Tracers >0 outside Greenland\n');
    fprintf('npnts = %i\n',length(I));
    error(' quitting ...');
  end
  
  szz = sprintf('min C=%5.2d, \nmax C=%5.2d',mnA,mxA);
  
  clf;
  contour(H,[0 0],'Color',[0.8 0.8 0.8]);
  hold on;
  pcolor(F); shading flat;
  caxis([0 100]);
  axis('equal');
  set(gca,'xlim',[250 450],'ylim',[250 550]);
  colorbar;
  spp=sprintf('Tracer conc., mo=%i, level=%i', im,kplt);
  title(spp);
  text(255,500,szz);
  drawnow

end;

% Plot relaxation scale:
fina = sprintf('%srelax_rmutr_Greenland.a',PTH.data);
finb = sprintf('%srelax_rmutr_Greenland.b',PTH.data);
fina_fid = fopen(fina,'r','ieee-be');
%finb_fid = fopen(finb,'wt');
dmm=fread(fina_fid,IJDM,'float32');  % read 2D field (1 layer)
%  dm1=fread(fid1,npad,'float32','ieee-be');  % read npad 

A=reshape(dmm,ID,JD);
A=A';

R=A(jg1:jg2,ig1:ig2);

figure(2); clf;
contour(H,[0 0],'Color',[.8 .8 .8]);
hold on;
R(R==0)=nan;
pcolor(R); shading flat;
%caxis([0 500]);
colorbar;
spp=sprintf('Relaxation scale, sec^{-1}');
title(spp);
drawnow

fclose(fina_fid);





