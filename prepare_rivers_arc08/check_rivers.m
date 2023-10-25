% Read rivers
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

%expt = 112; % NCAR + Greenland runoff
expt = 110; % river climatology
YY=2005;
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

if ~isempty(Iland),
  fprintf('ERR: Found %i river points on land\n',length(Iland));
  [jl,il]=ind2sub(size(HH),Iland);
  figure(2); clf;
  contour(HH,[0 0],'k');
  hold on;
  plot(IJarc(:,1),IJarc(:,2),'r.-');
  
  pcolor(tmp); shading flat;
  caxis([-10 -2]);
  colorbar
  
  plot(il,jl,'m*');
  title('River points on land');
  keyboard
end

f_gr=0; % plot greenland sources, km3/mo
if f_gr==1
  figure(10); clf;
  Lmsk = HH*0;
  Lmsk(HH<0)=1;
  cmL = [0.92 0.92 0.92; 1 1 1];

  pcolor(Lmsk); shading flat;
  colormap(cmL);
  freezeColors;

  %contour(HH,[0 0],'k');
  hold on;
  contour(HH,[-2000:250:-5],'Color',[0.6 0.6 0.6]);
  pcolor(Fm); shading flat;
  caxis([0 1.]);

  colormap(cmp);
  axis('equal');
  set(gca,'xlim',[550 630],...
	  'ylim',[675 775],...
	  'xtick',[],...
	  'ytick',[]);
  hc = colorbar;
  set(hc,'Fontsize',16,...
	 'TickLength', 0.025);
  tts=sprintf('Greenland FW forcing, 0.08 HYCOM, %i/%i',YY,MM);
  title(tts);
  
  btx='check_rivers.m';
  bottom_text(btx,'pwd',1);
  
end



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

Rarc = R-Rgr;
Rarc_ann = Rarc*3600*24*mdays*1e-9; % Arctic only, no Gr., km3/mo
Rgr_ann  = Rgr*3600*24*mdays*1e-9;  % Greenland, km3/mo
Rpc_ann  = Rpc*3600*24*mdays*1e-9;  % Pacific runoff, km3/mo
% Note that total runoff in Rgr_ann - monthly runoff
% should be similar to original Greenland runoff data
% compare with anls_GreenlandRunoff_v3.m - the only
% difference could be due to different versions
% of Bamber's data (2012 vs 2018)
% To compare: select year/month
% in anls_GreenlandRunoff_v3.m: FTMn(7,58) - gives 
% monthly total runoff for July year 2015
% FTMn(7,58) = 335.254 km3/mo
% compared to Rgr_ann for July 7 2015 = 335.465945887116
% good agreement


LMSK = HH*0;
LMSK(HH<0)=1;

figure(1);
clf
axes('Position',[0.1 0.1 0.85 0.8]);
pcolor(LMSK); shading flat;
cmp1=[0 0 0; 1 1 1];
colormap(cmp1);
freezeColors;
hold on;
contour(HH,[-1000 -1000],'Color',[0.7 0.7 0.7]);
axis('equal');
plot(IJarc(:,1),IJarc(:,2),'m.-');
plot(IJgr(:,1),IJgr(:,2),'g.-');
plot(IJpc(:,1),IJpc(:,2),'c.-');
pcolor(tmp); shading flat;
caxis([-8 -4]);
colormap('default');
colorbar
set(gca,'xlim',[10 nn],...
	'ylim',[10 mm-300]);
stt=sprintf('%s, log10(Runoff, m/s), M=%i',flriva(end-24:end),MM);
title(stt,'Interpreter','none');
txb = 'check_rivers.m';
bottom_text(txb,'pwd',1);

txt{1}=sprintf('Arctic(no Gr) Runoff, km3/mo: %7.1f',Rarc_ann);
txt{2}=sprintf('Greenl Runoff, km3/mo: %7.1f',Rgr_ann);
txt{3}=sprintf('N.Pacif Runoff, km3/mo: %7.1f',Rpc_ann);
axes('Position',[0.1 0.93 0.3 0.075]);
text(0.08,0.1,txt,'Fontsize',10);
set(gca,'xtick',[],'ytick',[],...
	'xlim',[0 1.5],'ylim',[0 0.25],...
	'box','off');

%title('log10(River runoff), m/s')

% Check if rivers are in the ocean points:
B2=HH;
J=find(B2==0);
B2(J)=NaN;
J2=find(~isnan(B2));
B2(J2)=0;
B2(J)=1;

J3=find(tmp~=0);
tmp(J3)=1;

D=tmp+B2;  % 
 I=find(D>1);
if ~isempty(I); error('Some river points are on land'); end


fclose(riv_fid);

%keyboard

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
if ~isempty(Iff), % no Greenland
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

  figure(3); clf;
  pcolor(Lmsk); shading flat;
  colormap(cmL);
  freezeColors;

  %contour(HH,[0 0],'k');
  hold on;
  contour(HH,[-5000:500:-5],'Color',[0.6 0.6 0.6]);
  %scatter(iff,jff,18,C,'filled');
  scatter(iff,jff,36,C,'filled');
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

end


% Plot selected regions for budget analysis
BX = sub_define_boxes(HH,LON,LAT,0);
for ib=1:5
  iBG = BX(ib).IJ;
  iBG(end+1,:)=iBG(1,:);
  plot(iBG(:,1),iBG(:,2),'r.-');
end


%pcolor(Fm); shading flat

ixx = max(strfind(flriva,'/'));
faa = flriva(ixx+1:end);
Ftot = nansum(nansum(Fm));
stl=sprintf('%s, FWFlux, km3/mo, IntgrFWF=%6.1f km3/mo %i/%2.2i',...
	    faa,Ftot,YY,MM);
title(stl,'Interpreter','none');
bottom_text(txb,'pwd',1);


IIk = find(A>0);
R = sum(A(IIk).*ACell(IIk))*3600*24*mdays*1e-9; %km3/mo
% No boxes - similar to ARCc0.04
figure(4); clf;
%pcolor(Lmsk); shading flat;
%colormap(cmL);
%freezeColors;

contour(HH,[0 0],'k');
hold on;
contour(HH,[-5000:500:-5],'Color',[0.6 0.6 0.6]);
%scatter(iff,jff,18,C,'filled');
scatter(iff,jff,36,C,'filled');
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
set(gca,'xlim',[520 975],...
	'ylim',[358 1098]);

set(gca,'xtick',[],...
	'ytick',[]);

stl=sprintf('%s, Runoff km3/mo, %i/%2.2i, Arc+Gr=%6.1f Greenl=%5.1f',...
	    faa,YY,MM,R,Ftot);
title(stl,'Interpreter','none');
bottom_text(txb,'pwd',1);


