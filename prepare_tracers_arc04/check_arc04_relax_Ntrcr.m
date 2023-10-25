% ARCc 0.04
% Tracers in hycom: prepared for N tracers
% Tracer concentration is calculated as
% m3/sec - river runoff(m/s) * grid cell area from ARCc0.08
% e.g., relax_tracer_Gr_AmEurRiv_Bering.m
%  Read in: for given month - read for tracer 1:ntr, all layers
% see forfun.f
% Check hycom binary prepared in relax_tracer_Gr_AmEurRiv_Bering.m
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close

regn   = 'ARCc0.04';
ntopo  = 17;   % HYCOM topo version
nmtopo = '17DD'; 
readall= 0; % read all records and check min/max
trplt  = 1;   % tracer to plot
kplt   = 1;    % layer to plot
mplt   = 5;    % month to plot
nTr    = trplt; 

% # of tracers in the file:
ntracr=5;  % # of tracers


fprintf('Checking tracers: %s\n',regn)
fprintf('HYCOM Topo=%2.2i, # Tracers=%i\n',ntopo,ntracr);
fprintf('Plotting: Tracer %i, Layer %i, Month %i\n\n',trplt,kplt,mplt);

PTH.data  = '/Net/mars/ddmitry/hycom/ARCc0.04/force/relax/';
PTH.topo  = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/topo_grid/';
PTH.river = '/Net/mars/ddmitry/hycom/ARCc0.04/force/rivers/';

switch(trplt)
 case(1)
  trnm='Greenland';
 case(2)
  trnm='MacKenzie';
 case(3)
  trnm='East EurRiv';
 case(4)
  trnm='West EurRiv';
 case(5)
  trnm='Bering';
end

% Get topo and grid:
fltopo=sprintf('%sdepth_ARCc0.04_%s.nc',PTH.topo,nmtopo);
fprintf('Reading topo: %s\n',fltopo);
HH   = nc_varget(fltopo,'Bathymetry');
LAT  = nc_varget(fltopo,'Latitude');
LON  = nc_varget(fltopo,'Longitude');
%[m,n]= size(HH);
[mm,nn]= size(HH);
[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2

% Greenland subdomain:
%jg1=385;
%jg2=1954;
%ig1=332;
%ig2=n;
jg1 = 1;
jg2 = mm;
ig1 = 1;
ig2 = nn;

%ftrca = sprintf('%srelax_trcr_Greenl_Riv_Bering.a',PTH.data);
%ftrcb = sprintf('%srelax_trcr_Greenl_Riv_Bering.b',PTH.data);
ftrca = sprintf('%srelax_trcr%2.2i_Greenl_Riv_Bering_T%2.2i.a',...
		PTH.data,ntracr,ntopo);
ftrcb = sprintf('%srelax_trcr%2.2i_Greenl_Riv_Bering_T%2.2i.b',...
		PTH.data,ntracr,ntopo);

fida = fopen(ftrca,'r');  % 
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

% Check all records, min/max
% Read all records
if readall>0
  for im=1:12
    fprintf('     MONTH: %i\n',im);
    for ktr=1:ntracr
      fprintf('     Tracer %i: \n',ktr);
      for k=1:nlev
	fprintf(':: Reading mo=%i, tracer=%i, lev=%i\n',im,ktr,k);
	dmm=fread(fida,IJDM,'float32','ieee-be');  % read 2D field (1 layer)
	dm1=fread(fida,npad,'float32','ieee-be');  % read npad 
	aa=fgetl(fidb);
	disp(aa);
	A=dmm;
	mnA=min(A);
	mnA0=0;
	I=find(A>0);
	if ~isempty(I)
	  mnA0=min(A(I));
	end
	mxA=max(A);
	fprintf('min C=%9.6f (>0: %9.6f), \nmax C=%9.6f\n',mnA,mnA0,mxA);
      end
    end
  end
end


ngrd=0;
%for kplt=1:41 % uncomment to read tracers in all layers
fprintf('Plotting mo=%i tracer=%i lev=%i \n',mplt,trplt,kplt);
nrec=(mplt-1)*(ntracr*nlev)+(trplt-1)*nlev+kplt-1; % this should match field order in dummy
stat=fseek(fida,nrec*(IJDM+npad)*4,-1);
if stat<0,
  error('Couldnot find record location in %s\n',ftrcaOLD);
end

dmm=fread(fida,IJDM,'float32','ieee-be');  % read 2D field (1 layer)
F=reshape(dmm,ID,JD)';
ll=find(F>0);
fprintf('# of river cells %i\n',length(ll));
%ngrd=ngrd+length(ll);
%fprintf('Total # of river cells  in all layers = %i\n',ngrd);
%end

%H=HH(jg1:jg2,ig1:ig2);
%F=A(jg1:jg2,ig1:ig2);

F(F==0)=nan;

mnA=min(min(F(F>0)));
mxA=max(max(F));

szz = sprintf('min C=%5.2d, \nmax C=%5.2d',mnA,mxA);

% Plotting tracer - should be same as in ARCc0.08
% 
%c1=0;
%c2=50;
%CMP = create_colormap2_1(100,c1,c2);
%cmp = CMP.colormap;
%cnt = CMP.intervals;

TrMn = 1e-10;
F(F<=TrMn) = nan;
lTr = log(F);

% Total Mass in 1 layer, GT
M = nansum(nansum(F.*Acell))*1e-12;

nf = 1;
xlim1 = 0;
xlim2 = nn;
ylim1 = 0;
ylim2 = mm;
stl=sprintf('%s, T=%s, M=%6.3f GT, Tr %i (%s), mo=%i, level=%i', ...
	    regn,nmtopo,M,trplt,trnm,mplt,kplt);
sub_plot_tracers(lTr,nf,HH,xlim1,xlim2,...
		 ylim1,ylim2,LON,LAT,stl,nTr,...
		 'c1',-3,'c2',2);

%title(spp,'fontsize',12);
text(800,100,szz);
txb = 'check_arc04_relax_Ntrcr.m';
bottom_text(txb,'pwd',1);

fclose(fida);
fclose(fidb);


% Plot relaxation time scale:
ftrca = sprintf('%srelax_rmutr_Greenl_Riv_Bering_T%2.2i.a',PTH.data,ntopo);
ftrcb = sprintf('%srelax_rmutr_Greenl_Riv_Bering_T%2.2i.b',PTH.data,ntopo);

fida = fopen(ftrca,'r');  % 
dmm=fread(fida,IJDM,'float32','ieee-be');  % read 2D field (1 layer)
RL=reshape(dmm,ID,JD)';

RL(RL==0)=nan;
%RL=A(jg1:jg2,ig1:ig2);
RL=1./RL*(1/3600); % sec -> hours

c1=0;
c2=24;
CMP = create_colormap2_3(100,c1,c2);
cmp = CMP.colormap;
cnt = CMP.intervals;


keyboard

figure(2);
clf;
contour(HH,[0 0],'Color',[0.8 0.8 0.8]);
hold on;
pcolor(RL); shading flat;
%  caxis([0  c2]);
axis('equal');
%set(gca,'xlim',[400 nn],'ylim',[800 4400]);
caxis([c1 c2]);
colormap(cmp);

colorbar;
Is=max(strfind(ftrca,'/'));
spp=sprintf('%s, Relax Time (hrs) %s',regn,ftrca(52:end));
title(spp,'interpreter','none','Fontsize',12);

txb = 'check_relax_Ntrcr.m';
bottom_text(txb,'pwd',1);
drawnow

fclose(fida);






