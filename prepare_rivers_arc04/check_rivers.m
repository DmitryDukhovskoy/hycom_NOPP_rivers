% Read rivers
% Units:  m/s
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom
startup;

close all
clear

TV = '17DD';  % topo name
ntopo = 17;
PTH.data = '/Net/mars/ddmitry/hycom/ARCc0.04/force/rivers/';
PTH.topo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/topo_grid/';
PTH.dataout = '/Net/mars/ddmitry/hycom/ARCc0.04/force/riversNCAR/';
%flt=[pth,'regional.grid.b'];

YearGr = 2016;
YY=YearGr;
MM     = 7;  % month to check
YYx=YearGr;
MMx=MM+1;
if MMx==13,
  MMx=1;
  YYx=YY+1;
end
mdays = datenum(YYx,MMx,1)-datenum(YY,MM,1);

%flriva = sprintf('%srivers_%s.a',PTH.data,TV); % river clim, no Greenland
%flrivb = sprintf('%srivers_%s.b',PTH.data,TV);
%flriva = sprintf('%srivers_%s_Greenland_%4.4i.a',PTH.data,TV,YearGr); 
%flrivb = sprintf('%srivers_%s_Greenland_%4.4i.b',PTH.data,TV,YearGr);
flriva=sprintf('%srivers_%2.2i_NCAR_Gr_%4.4i.a',PTH.dataout,ntopo,YearGr); 
flrivb=sprintf('%srivers_%2.2i_NCAR_Gr_%4.4i.b',PTH.dataout,ntopo,YearGr);
fltopo = sprintf('%sdepth_ARCc0.04_%s.nc',PTH.topo,TV);

fprintf('Reading %s\n',flriva);
fprintf('Year %i, Mo=%i\n',YearGr,MM);

% Get topo and grid:
HH   = nc_varget(fltopo,'Bathymetry');
LAT = nc_varget(fltopo,'Latitude');
LON = nc_varget(fltopo,'Longitude');
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
IJarc=IJarc*2;
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
IJgr=IJgr*2;
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
IJpc=IJpc*2;
dmm = inpolygon(II,JJ,IJpc(:,1),IJpc(:,2));
INpc = find(dmm==1);
OUTpc = find(dmm==0);


% Reading rivers & bathymetry:
IDM=nn;
JDM=mm;
IJDM=IDM*JDM;
npad=4096-mod(IJDM,4096);

% Rivers:
fprintf('Reading rivers %s\n',flriva);
riv_fid=fopen(flriva,'r');
fseek(riv_fid,(MM-1)*4*(npad+IJDM),-1);
[A,counta]=fread(riv_fid,IJDM,'float32','ieee-be');
  toto = fread(riv_fid,npad,'float32','ieee-be');
%  dmm=fread(riv_fid,1,'int','ieee-be'); <-- Do not need this, direct access format
%I=find(A==0);
%A(I)=NaN;
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
  
  pcolor(tmp); shading flat;
  caxis([-10 -2]);
  colorbar
  
  plot(il,jl,'m*');
  title('River points on land');
  keyboard
end


% Check total river runoff, m3/s, m3/mo
IIk = find(A>0);
%R = sum(A(IIk).*ACell(IIk));
R = sum(A(IIk).*ACell(IIk))*3600*24*mdays*1e-9; %km3/mo
RGr = nansum(nansum(Fm));

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


% Plot Greenlnd
xl1=1040;
xl2=1950;
yl1=715;
yl2=2195;


figure(3);
clf
contour(HH,[0 0],'Color',[0. 0. 0.]);
hold on;
%contour(HH,[-100:20:-10],'Color',[0.9 0.9 0.9]);
contour(HH,[-5000:500:-10],'Color',[0.6 0.6 0.6]);
%pcolor(tmp); shading flat;
%pcolor(Fm); shading flat;
scatter(iff,jff,36,C,'filled');
colormap(cmp);
caxis([c1 c2]);

axis('equal');
set(gca,'xlim',[xl1 xl2],...
	'ylim',[yl1 yl2],...
	  'xtick',[],...
	  'ytick',[]);

hc = colorbar;
set(hc,'Fontsize',16,...
       'TickLength', 0.025);

LL = LON;
LL(2200:end,:)=nan;
contour(LL,[-180:10:180],'Color',[0.6 0.6 0.6]);
contour(LAT,[30:10:89],'Color',[0.6 0.6 0.6]);


stt=sprintf('ARCc0.04, T=%s, Runoff, km3/mo, %i/%i, Arc+Gr=%6.1f, Greenl=%5.1f',...
	    TV,YearGr,MM,R,RGr);
title(stt);
btx = 'check_rivers.m';
bottom_text(btx,'pwd',1);

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







