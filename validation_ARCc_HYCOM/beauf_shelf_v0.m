% Check T/S water mass profiles/ water masses
% from AO 0.08 HYCOM-CICE
% with improved rivers and Greenland runoff on
%
% Beaufort shelf
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

shlf='Beaufort';
expt = 112;  % epxeriment with Greenland runoff and monthly Arctic rivers
regn = 'ARCc0.08';


rg=9806;  % convert pressure to depth, m
hg     = 2^100;  % "huge" in HYCOM used for land masking ONLY!
huge   = hg;     % 0-depth values are not = huge
hgg=1e20; 
btx = 'beauf_shelf.m';


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

% define Beauf. shelf in ARCc0.08 grid:
IJ=[394        1552
         449        1507
         437        1540
         445        1574
         442        1607
         456        1630
         507        1665
         518        1675
         534        1680
         548        1690
         645        1687
         612        1756
         527        1724
         416        1659
         395        1640
         377        1569];



% Read obs:
fobs=sprintf('%stsz_WIDE_20558.dat',pthdata);
fprintf('Reading obs %s\n',fobs);
fidb=fopen(fobs,'r');
frewind(fidb);
%aa=fgetl(fidb);
S=textscan(fidb,'%u%f%f%f%f%f%f');
fclose(fidb);

% Group by observations
cc=0;
ID=S{1};
TM=S{2}+datenum(1900,1,1)-1;
LTo=S{3};
LNo=S{4};
ZZ=S{5};
TT=S{6};
SS=S{7};
TT(TT>25)=nan;
TT(TT<-1.89)=nan;
SS(SS>36)=nan;
SS(SS<1)=nan;
I=find(LNo>180);
LNo(I)=LNo(I)-360;
Iz=find(ZZ<2); % surface data
SS(Iz)=nan;
TT(Iz)=nan;

clear OBS DV
nk=length(ID);
iend=0;
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
  
% find closest HYCOM index for LOT/LAT:
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
save(fobsmat,'II','JJ');


figure(1); clf;
contour(HH,[0 0],'k');
hold
contour(HH,[-5000:1000:-900],'Color',[0.7 0.7 0.7]);
contour(HH,[-500:100:-90],'Color',[0.85 0.85 0.85]);
axis('equal');
set(gca,'xlim',[350 700],...
	'ylim',[1400 1800]);
plot(IJ(:,1),IJ(:,2),'r-');


plot(II,JJ,'b.');
stl=sprintf('Obs. Beaufort Shelf, Summer, %i-%i',min(DV(:,1)),max(DV(:,1)));
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


Zbin=[-1000; -100; -50; -30; 0];
nb=length(Zbin);
CLR=colormap(parula(nb-1));
for iz=1:nb
  lbb{iz}=sprintf('%i',Zbin(iz));
end

dx=0.4;
POS=[0.08 0.52 dx dx; ...
     0.55 0.52 dx dx; ...
     0.08 0.08 dx dx; ...
     0.55 0.08 dx dx];

% Plot by depth ranges:
figure(2); clf;
axes('Position',[0.08 0.42 0.74 0.5]);
hold on;
for ik=1:nobs
  SS=OBS(ik).S;
  TT=OBS(ik).T;
  ZZ=OBS(ik).Z;
  for iz=1:nb-1
    iz1=Zbin(iz);
    iz2=Zbin(iz+1);
    I=find(ZZ>=iz1 & ZZ<iz2);
    if isempty(I); continue; end;
    clr=CLR(iz,:);
    plot(SS(I),TT(I),'.','Color',clr);
  end
  
end
contour(SGR,TGR,sgm0,[0:2:28],'k','linewidth',1.5);

colormap(CLR);
caxis([1 nb]);
cp=colorbar;
set(cp,'Position',[0.85 0.42 0.02 0.5],...
       'Ticks',Zbin,...
       'TickLength',0.02,...
       'Ticks',[1:nb],...
       'Ticklabels',lbb)


axis('equal');
set(gca,'tickdir','out',...
	'xlim',[24 35.5],...
	'xtick',[0:2:40],...
	'ylim',[-2 12],...
	'ytick',[-2:2:16],...
	'Fontsize',14);
xlabel('S');
ylabel('T');

stl=sprintf('Obs. Beaufort Shelf, Summer, %i-%i',min(DV(:,1)),max(DV(:,1)));
title(stl);

bottom_text(btx,'pwd',1,'Position',[0.08 0.2 0.4 0.04]);


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
  ZM=[ZM;SHLF.ZM];
  Th=[Th;SHLF.Temp];
  Sh=[Sh;SHLF.Salin];
end

%I=find(ZM>=-5);
%ZM(I)=10;

figure(3); clf;
axes('Position',[0.08 0.42 0.74 0.5]);
hold on;
for iz=nb-1:-1:1
  iz1=Zbin(iz);
  iz2=Zbin(iz+1);
  I=find(ZM>=iz1 & ZM<iz2);
  if isempty(I); continue; end;
  clr=CLR(iz,:);
  plot(Sh(I),Th(I),'.','Color',clr);
end
  
contour(SGR,TGR,sgm0,[0:2:28],'k','linewidth',1.5);

colormap(CLR);
caxis([1 nb]);
cp=colorbar;
set(cp,'Position',[0.85 0.42 0.02 0.5],...
       'Ticks',Zbin,...
       'TickLength',0.02,...
       'Ticks',[1:nb],...
       'Ticklabels',lbb)


axis('equal');
set(gca,'tickdir','out',...
	'xlim',[24 35.5],...
	'xtick',[0:2:40],...
	'ylim',[-2 12],...
	'ytick',[-2:2:16],...
	'Fontsize',14);
xlabel('S');
ylabel('T');

stl=sprintf('HYCOM 0.08-112 Beaufort Shelf, Summer, %i-%i',min(DV(:,1)),max(DV(:,1)));
title(stl);

bottom_text(btx,'pwd',1,'Position',[0.08 0.2 0.4 0.04]);
