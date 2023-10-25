% Plot regino map
% with bath + Gr runoff
% + obs points
% Plot river sources
% Check onland river sources
% Units:  m/s
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
startup;

close all
clear

PTH.data1='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/force/rivers/';
PTH.data = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/force/riversNCAR/';
PTH.topo='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
%flt=[pth,'regional.grid.b'];

expt = 112; % NCAR + Greenland runoff
%expt = 110; % river climatology
YY=2015;
MM = 7; % month
YYx=YY;
MMx=MM+1;
if MMx==13,
  MMx=1;
  YYx=YY+1;
end

mdays = datenum(YYx,MMx,1)-datenum(YY,MM,1);

%flriva=sprintf('%srivers_09.a',PTH.data);
%flrivb=sprintf('%srivers_09.b',PTH.data);
%flriva=sprintf('%srivers_11_Greenland_%i.a',PTH.data,YY);
%flrivb=sprintf('%srivers_11_Greenland_%i.b',PTH.data,YY);
switch(expt)
 case(110);
  flriva=sprintf('%srivers_11.a',PTH.data1);  % expt 11.0 clim. rivers, no Gr.
  flrivb=sprintf('%srivers_11.b',PTH.data1);
 case(112)
  flriva=sprintf('%srivers_11_NCAR_Gr_%4.4i.a',PTH.data,YY); % NCAR riv+Green.
  flrivb=sprintf('%srivers_11_NCAR_Gr_%4.4i.b',PTH.data,YY); % NCAR riv+Green.
end

%fltopo=sprintf('%sdepth_ARCc0.08_09.nc',PTH.topo);
fltopo=sprintf('%sdepth_ARCc0.08_11.nc',PTH.topo);

fprintf('Reading %s\n',flriva);
fprintf('Expt=%3.3i, Year %i, Mo=%i\n',expt,YY,MM);

% Get topo and grid:
HH   = nc_varget(fltopo,'Bathymetry');
LAT = nc_varget(fltopo,'Latitude');
LON = nc_varget(fltopo,'Longitude');
[m,n]= size(HH);
[mm,nn]= size(HH);
[DX,DY]=sub_dx_dy(LON,LAT);
ACell=DX.*DY;

% Arctic Domain with Greenland
IJarc=[  380         353
         505         355
	 705         376
	 841         558
	 937         575
        1184         449
        1298         686
        1483         612
        1594         664
        1594        1323
	1235        1916
         443        1916
         188        1398
          51        1005
          44         531
         147         234];
IJarc(end+1,:)=IJarc(1,:);
[II,JJ]=meshgrid((1:nn),(1:mm));
dmm = inpolygon(II,JJ,IJarc(:,1),IJarc(:,2));
IN = find(dmm==1);
OUT = find(dmm==0);

% Greenland: 
IJgr =[  602         372
         912         645
         982         960
         917        1080
         785        1076
         583        1084
         558         923
         503         488
         539         401];
IJgr(end+1,:)=IJgr(1,:);
dmm = inpolygon(II,JJ,IJgr(:,1),IJgr(:,2));
INg = find(dmm==1);
OUTg = find(dmm==0);

% Pacific:
IJpc = [           6        2492
          30        1930
         331        1870
         631        1926
        1141        1948
        1325        2235
        1415        2500];
IJpc(end+1,:)=IJpc(1,:);
dmm = inpolygon(II,JJ,IJpc(:,1),IJpc(:,2));
INpc = find(dmm==1);
OUTpc = find(dmm==0);


% Reading rivers & bathymetry:
IDM=n;
JDM=m;
IJDM=IDM*JDM;
npad=4096-mod(IJDM,4096);

% Rivers:
riv_fid=fopen(flriva,'r');
fseek(riv_fid,(MM-1)*4*(npad+IJDM),-1);
[A,counta]=fread(riv_fid,IJDM,'float32','ieee-be');
toto = fread(riv_fid,npad,'float32','ieee-be');
A=reshape(A,IDM,JDM)';

[jj,ii]=find(A>0);

% Check river points on land:
Iland=find(A>0 & HH>=0);

Fm = A.*ACell*3600*24*mdays*1e-9; %km3/mo
Fm(A<1e-23)=nan;
Fm(OUTg) = nan;

tmp=log10(A);
tmp(A==0)=NaN;


% Check total river runoff, m3/s, m3/mo
% Arctic Ocean+ Greenland:
dmm = A;
dmm(OUT) = nan;
II = find(dmm>0);
R = sum(A(II).*ACell(II)); % m/s -> m3/s
dmm = A;
dmm(OUTg) = nan;
II = find(dmm>0);
Rgr = sum(A(II).*ACell(II)); % m/s -> m3/s
dmm = A;
dmm(OUTpc) = nan;
II = find(dmm>0);
Rpc = sum(A(II).*ACell(II)); % m/s -> m3/s
clear dmm


% =========================
% Map
% Plot individual FW sources in Greenland
% =========================
IGR=[    608        1042
         543         941
         520         652
         523         454
         602         400
         909         653
         955         976
         896        1069
         768        1070];

[X,Y] = meshgrid([1:mm],[1:nn]);
%IBB = inpolygon(X,Y,IGR(:,1),IGR(:,2));
%Iout = find(IBB==0);
%Iin  = find(IBB==1);

%Fm(Iout)=nan;

c1=0;
c2=1.4;
ncc = 200;
cmp=colormap(parula(ncc));
cnt = (c1:(c2-c1)/ncc:c2);


Iff=find(~isnan(Fm));
[jff,iff]=ind2sub(size(HH),Iff);
nss = length(Iff);
clear C
for ik=1:nss
  i0 = Iff(ik);
  fvl = Fm(i0);
  ix = max(find(cnt<=fvl));
  ix=min([ncc,ix]);
  C(ik,:) = cmp(ix,:);
end

Lmsk = HH*0;
Lmsk(HH<0)=1;
cmL = [0.4 0.4 0.4; 1 1 1];

% Get observational sites from I. Yashayev
POBS = sub_obs_locations;



figure(3); clf;
pcolor(Lmsk); shading flat;
colormap(cmL);
freezeColors;

%contour(HH,[0 0],'k');
hold on;
contour(HH,[-5000:500:-5],'Color',[0.6 0.6 0.6]);
%scatter(iff,jff,18,C,'filled');
scatter(iff,jff,86,C,'filled');
caxis([c1 c2]);

colormap(cmp);
hc = colorbar;
set(hc,'Fontsize',16,...
       'TickLength', 0.025);
  
LL = LON;
LL(1100:end,:)=nan;
contour(LL,[-180:10:180],'Color',[0.6 0.6 0.6]);
contour(LAT,[30:10:89],'Color',[0.6 0.6 0.6]);

axis('equal');
%set(gca,'xlim',[300 1250],...
%	'ylim',[150 1100]);
set(gca,'xlim',[200 1350],...
	'ylim',[1 1100]);

set(gca,'xtick',[],...
	'ytick',[]);


% Plot locations of in situ obs - Yashayaev
nobs = length(POBS);
for ii=1:nobs
  X=POBS(ii).Lon;
  Y=POBS(ii).Lat;
  X=X(:);
  Y=Y(:);
  IJ = sub_indices_domain([X,Y],LON,LAT); 
  if length(X)>1,
    plot(IJ(:,1),IJ(:,2),'r-','Color',[0.6 0 0],'Linewidth',1.6);
  else
    plot(IJ(1),IJ(2),'r.','Color',[0.6 0 0],'Markersize',18);
  end
  
end  


% Plot my selected regions for budget analysis
%BX = sub_define_boxes(HH,LON,LAT,0);
%for ib=1:5
%  iBG = BX(ib).IJ;
%  iBG(end+1,:)=iBG(1,:);
%  plot(iBG(:,1),iBG(:,2),'r.-');
%end


%pcolor(Fm); shading flat
btx = 'Greenl_bath_runoff_maps.m';
ixx = max(strfind(flriva,'/'));
faa = flriva(ixx+1:end);
Ftot = nansum(nansum(Fm));
stl=sprintf('%s, FWFlux, km3/mo, IntgrFWF=%6.1f km3/mo %i/%2.2i',...
	    faa,Ftot,YY,MM);
title(stl,'Interpreter','none');
bottom_text(btx,'pwd',1,'Position',[0.02 0.1 0.4 0.05]);



