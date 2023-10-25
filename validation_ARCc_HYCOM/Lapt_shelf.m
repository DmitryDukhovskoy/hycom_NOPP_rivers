% Check T/S water mass profiles/ water masses
% from AO 0.08 HYCOM-CICE
% with improved rivers and Greenland runoff on
%
% Lapttev shelf
% Data downloaded from
% http://oregon.iarc.uaf.edu/db/tszresultR.cgi
%
% File format:
% WIDE File contains no header information or comments, 
% and each row is tab-delimited containing:
% INDEX TIME LONG(X) LAT(Y) DEPTH TEMPERATURE SALINITY
% where "TIME" counts the number of julian days since 1900-01-01.
%
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
%addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_gom04/model_nest_ob;
startup

shlf='Lapt';
expt = 112;  % epxeriment with Greenland runoff and monthly Arctic rivers
regn = 'ARCc0.08';


rg=9806;  % convert pressure to depth, m
hg     = 2^100;  % "huge" in HYCOM used for land masking ONLY!
huge   = hg;     % 0-depth values are not = huge
hgg=1e20; 
btx = 'Lapt_shelf.m';


pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthdata = '/Net/data2/ddmitry/Arctic_Data/BeaufortShelf/';
pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat2/',expt);

% -------------------------
% Get grid and bath:
% My bathymetry, Model bathymetry:
% -------------------------
ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

% define  shelf in ARCc0.08 grid:
IJ=[1120        1571
        1130        1641
        1147        1708
        1179        1756
        1216        1769
        1206        1696
        1266        1664
        1274        1557
        1247        1541
        1197        1517
        1182        1565
        1152        1581
        1134        1571];


% Read 3 files with obs:
fobs=sprintf('%sLaptSea.dat',pthdata);
fprintf('Reading obs %s\n',fobs);
fidb=fopen(fobs,'r');
frewind(fidb);
%aa=fgetl(fidb);
S=textscan(fidb,'%f%f%f%f%f%f%f');
fclose(fidb);

% Group by observations
ID=S{1};
TM=S{2}+datenum(1900,1,1)-1;
LTo=S{3};
LNo=S{4};
ZZ=S{5};
TT=S{6};
SS=S{7};

fobs=sprintf('%sLaptSea2.dat',pthdata);
fprintf('Reading obs %s\n',fobs);
fidb=fopen(fobs,'r');
frewind(fidb);
%aa=fgetl(fidb);
S=textscan(fidb,'%f%f%f%f%f%f%f');
fclose(fidb);

ID=[ID;S{1}];
dmm=S{2}+datenum(1900,1,1)-1;
TM=[TM;dmm];
LTo=[LTo;S{3}];
LNo=[LNo;S{4}];
ZZ=[ZZ;S{5}];
TT=[TT;S{6}];
SS=[SS;S{7}];

fobs=sprintf('%sLaptSea3.dat',pthdata);
fprintf('Reading obs %s\n',fobs);
fidb=fopen(fobs,'r');
frewind(fidb);
%aa=fgetl(fidb);
S=textscan(fidb,'%f%f%f%f%f%f%f');
fclose(fidb);

ID=[ID;S{1}];
dmm=S{2}+datenum(1900,1,1)-1;
TM=[TM;dmm];
LTo=[LTo;S{3}];
LNo=[LNo;S{4}];
ZZ=[ZZ;S{5}];
TT=[TT;S{6}];
SS=[SS;S{7}];

fobs=sprintf('%sLaptSea4.dat',pthdata);
fprintf('Reading obs %s\n',fobs);
fidb=fopen(fobs,'r');
frewind(fidb);
%aa=fgetl(fidb);
S=textscan(fidb,'%f%f%f%f%f%f%f');
fclose(fidb);

% check overlap
aa=S{1};
for k=1:length(aa);
  ib=aa(k);
  J=find(ID==ib);
  if ~isempty(J)
    aa(k)=nan;
  end
end
Jn=find(~isnan(aa));
ID=[ID;aa(Jn)];
dmm=S{2}+datenum(1900,1,1)-1;
TM=[TM;dmm(Jn)];
aa=S{3};
LTo=[LTo;aa(Jn)];
aa=S{4};
LNo=[LNo;aa(Jn)];
aa=S{5};
ZZ=[ZZ;aa(Jn)];
aa=S{6};
TT=[TT;aa(Jn)];
aa=S{7};
SS=[SS;aa(Jn)];



TT(TT>25)=nan;
TT(TT<-1.89)=nan;
SS(SS>36)=nan;
SS(SS<1)=nan;
I=find(LNo>180);
LNo(I)=LNo(I)-360;
%Iz=find(ZZ<2); % surface data
%SS(Iz)=nan;
%TT(Iz)=nan;

clear OBS DV
nk=length(ID);
iend=0;
cc=0;
while iend<nk
  ist=iend+1;
  id=ID(ist);
  I=find(ID==id);
  iend=max(I);
  
% only observations on the shelf:  
  ln0=LNo(ist);
  lt0=LTo(ist);
% find closest HYCOM index for LOT/LAT:
  [ii,jj]=sub_find_indxHYCOM(LON,LAT,ln0,lt0);
  IP=inpolygon(ii,jj,IJ(:,1),IJ(:,2));
  
  if ~IP, continue; end;
  
  cc=cc+1;
  fprintf('cc=%i, Station: %i, %4.2f, %4.2f\n',cc, id,ln0,lt0);
  OBS(cc).ID=id;
  OBS(cc).TM=TM(ist);
  OBS(cc).LON=ln0;
  OBS(cc).LAT=lt0;
  OBS(cc).Z=-ZZ(ist:iend);
  OBS(cc).T=TT(ist:iend);
  OBS(cc).S=SS(ist:iend);
  dv=datevec(TM(cc));
  DV(cc,1:5)=dv(1:5);
  
% find closest HYCOM index for LON/LAT:
%  [ii,jj]=sub_find_indxHYCOM(LON,LAT,ln0,lt0);
  OBS(cc).HYCOM_I=ii;
  OBS(cc).HYCOM_J=jj;
  
end

nobs=cc;
% Plot stations:
clear II JJ
for ik=1:nobs
  II(ik)=OBS(ik).HYCOM_I;
  JJ(ik)=OBS(ik).HYCOM_J;
end

% Save obs locations:
% For extracting hycom at these locations
yr1=min(DV(:,1));
yr2=max(DV(:,1));
fobsmat=sprintf('%s%s_locations_%i-%i.mat',pthmat,shlf,yr1,yr2);
fprintf('Saving Obs. locations: %s\n',fobsmat);
save(fobsmat,'II','JJ','IJ');


figure(1); clf;
contour(HH,[0 0],'k');
hold
contour(HH,[-5000:1000:-900],'Color',[0.7 0.7 0.7]);
contour(HH,[-500:100:-90],'Color',[0.85 0.85 0.85]);
axis('equal');
set(gca,'xlim',[950 1310],...
	'ylim',[1400 1800]);
plot(IJ(:,1),IJ(:,2),'r-');
plot(II,JJ,'b.');
stl=sprintf('Obs. Laptev Shelf, Summer, %i-%i',min(DV(:,1)),max(DV(:,1)));
title(stl);


% Select Winter Water masses:
Iw=find(DV(:,2)>=10 | DV(:,2)<=3); % no winter obs.
Is=find(DV(:,2)>=6 & DV(:,2)<=9);

% TS diagrams
% Plot sigmas:
SGR=[0:0.1:35];
TGR=[-1.8:0.1:16];
ns=length(SGR);
nt=length(TGR);
clear sgm0
for ik=1:nt;  % 
  rho=sw_dens0(SGR,ones(1,ns)*TGR(ik))-1000;
  sgm0(ik,1:length(rho))=rho;
end;


%Zbin=[0; -30; -50; -100; -1000];
Zbin=[0; -30; -100; -1000];
nb=length(Zbin);
CLR=colormap(parula(nb-1));
for iz=1:nb
  lbb{iz}=sprintf('%i',Zbin(iz));
end

dx=0.38;
POS=[0.08 0.56 dx dx; ...
     0.55 0.56 dx dx; ...
     0.08 0.08 dx dx; ...
     0.55 0.08 dx dx];

% Show near-coast obs, 
% due to presence of ice in the model
% HYCOM does not have warm saline waters in the offshore 
% shelf region
% Plot by depth ranges:
figure(2); clf;
for iz=1:2
  xl1=24;
  xl2=35.5;
  dxx=2;
  if iz==1,
    xl1=10;
    dxx=5;
  end
  pos=POS(iz,:);
  axes('Position',pos);
  hold on;
  
  iz1=Zbin(iz+1);
  iz2=Zbin(iz);
  
  for ik=1:nobs
    SS=OBS(ik).S;
    TT=OBS(ik).T;
    ZZ=OBS(ik).Z;
    lat=OBS(ik).LAT;
    I=find(ZZ>=iz1 & ZZ<iz2);
    if iz2>-30 & lat>73.5; continue; end; % select near-coastal points only for shallow locations
    if isempty(I); continue; end;
    plot(SS(I),TT(I),'.','Color',[0.5 0.5 0.5]);
  end
  [c,h]=contour(SGR,TGR,sgm0,[0:2:28],'k','linewidth',1.5);
  clabel(c,h,'Fontsize',14,'LabelSpacing',400,'FontWeight','bold');
 
%  axis('equal');
  set(gca,'tickdir','out',...
	'xlim',[xl1 xl2],...
	'xtick',[0:dxx:40],...
	'ylim',[-2 12],...
	'ytick',[-2:2:16],...
	'Fontsize',14);
  xlabel('S');
  ylabel('T');

  stl=sprintf('Obs.Shlf %s, %4.1f-%4.1fm, Summer, %i-%i',...
	      shlf,abs(iz1),abs(iz2),min(DV(:,1)),max(DV(:,1)));
  title(stl,'Fontsize',14);
  
end

bottom_text(btx,'pwd',1);

%keyboard
  

%
% Plot hycom data:
% Subsampled into obs locations:
TMh=[];
ZM=[];
Th=[];
Sh=[];
for yr=yr1:yr2
%  fmat=sprintf('%s%s-%3.3i_TSprf_shelf_%s_%i.mat',pthmat,regn,expt,shlf,yr);
  fmat=sprintf('%s%s-%3.3i_TSprf_obsloc_%s_%i.mat',pthmat,regn,expt,shlf,yr);
  fprintf('Loading %s\n',fmat);
  load(fmat);
  
  TMh=[TMh;SHLF.TM'];  
  Ish=SHLF.IJ_shelf;
% Obs locatiiions  
%  ZM=[ZM;SHLF.ZM];
%  Th=[Th;SHLF.Temp];
%  Sh=[Sh;SHLF.Salin];
% 
% Whole shelf:
  ZM=[ZM;SHLF.ZM_sh];
  Th=[Th;SHLF.Temp_sh];
  Sh=[Sh;SHLF.Salin_sh];
end

%I=find(ZM>=-5);
%ZM(I)=10;
LATsh=LAT(Ish);

figure(3); clf;
for iz=1:2
  xl1=24;
  xl2=35.5;
  dxx=2;
  if iz==1,
    xl1=10;
    dxx=5;
  end
  
  pos=POS(iz,:);
  axes('Position',pos);
  hold on;
  
  iz1=Zbin(iz+1);
  iz2=Zbin(iz);
  
  I=find(ZM>=iz1 & ZM<iz2);
  ss=Sh(I);
  tt=Th(I);
  if iz1<-35
    J=find(LATsh<77);
    dmm=Sh;
    dmm(:,:,J)=nan;
    ss=dmm(I);
    dmm=Th;
    dmm(:,:,J)=nan;
    tt=dmm(I);
  end
  
  if isempty(I); continue; end;
  plot(ss,tt,'.','Color',[0.5 0.5 0.5]);

  [c,h]=contour(SGR,TGR,sgm0,[0:2:28],'k','linewidth',1.5);
  clabel(c,h,'Fontsize',14,'LabelSpacing',400,'FontWeight','bold');
 
%  axis('equal');
  set(gca,'tickdir','out',...
	'xlim',[xl1 xl2],...
	'xtick',[0:dxx:40],...
	'ylim',[-2 12],...
	'ytick',[-2:2:16],...
	'Fontsize',14);
  xlabel('S');
  ylabel('T');

  stl=sprintf('HYCOM 008-012 Shlf %s, %4.1f-%4.1fm, Summer, %i-%i',...
	      shlf,abs(iz1),abs(iz2),min(DV(:,1)),max(DV(:,1)));
  title(stl,'Fontsize',12);
end

bottom_text(btx,'pwd',1);
