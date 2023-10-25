% Monthly mean S from
% Aquarius
% Data provided by Subra and Rachel
%
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

% select a satellite:
%satl='Aquarius'; % 2012-2015
%satl='SMOS';
satl='SMAP';  % 2015-2018
%satl='SMOS';


iyr1 = 2015;
iyr2 = 2016;
%iyr2 = 2014;
%mplt = [12,1,2]; % average over N months, if not, mplt=1month
%mplt = [6,7,8]; % average over N months, if not, mplt=1month
mplt = [6,7,8]; % average over N months, if not, mplt=1month
nav = length(mplt);

pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = '/Net/gleam/dmitry/Aquarius/';

fmat = sprintf('%sDmitry_%s.mat',pthmat,satl);

fprintf('Loading %s\n',fmat);
AQ=load(fmat);

%ftopo = sprintf('%s/depth_%s_17DD.nc',pthtopo,regn); % 
ftopo = sprintf('%s/depth_ARCc0.08_11.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

SM={'J','F','M','A','M','J','J','A','S','O','N','D'};

TM=AQ.time;
DV=datevec(TM);
%S=AQ.sal;
slat=AQ.lat;
slon=AQ.lon;
m2=length(slat);
n2=length(slon);
nmm=m2*n2;

%ss=squeeze(S(:,:,10));
%ss=squeeze(S(:,:,11));
%ss=squeeze(S(:,:,12));
%ss=squeeze(S(:,:,13));

% Interpolate into HYCOM,
% subsample region to expedite calcu
bj1=1;
bj2=1300;
bi1=350;
bi2=1050;
sblon=LON(bj1:bj2,bi1:bi2);
sblat=LAT(bj1:bj2,bi1:bi2);

% HYCOM indices for plotting
findx=sprintf('%s%s_008hycom_indx.mat',pthmat,satl);
kpp=0;
if kpp==1
  for jj=1:m2
    for ii=1:n2
      kpp=kpp+1;
      if mod(kpp,1000)==0,
	fprintf('===> %4.1f%%\n',kpp/nmm*100);
      end

      x0=slon(ii);
      y0=slat(jj);
      dst=distance_spheric_coord(sblat,sblon,y0,x0);;
      [j1,i1]=find(dst==min(min(dst)),1);
      I_hycom(jj,ii)=i1+bi1-1;
      J_hycom(jj,ii)=j1+bj1-1;
    end
  end
  fprintf('Saving HYCOM indx %s\n',findx);
  save(findx,'I_hycom','J_hycom');
else
  fprintf('Loading HYCOM indx %s\n',findx);
  load(findx);
end


% Find summer mean S
smm=zeros(m2,n2);
ncc=smm;
for iyy=iyr1:iyr2
  for imm=1:nav
    im=mplt(imm);
    dnmb=datenum(iyy,im,15);
    itime=find(TM==dnmb);
    if isempty(itime),
      fprintf('No record: %s\n',datestr(TM(itime)));
      continue
    end
    if isfield(AQ,'sal')
      A=squeeze(AQ.sal(:,:,itime));
    elseif isfield(AQ,'SMAP')
      A=squeeze(AQ.SMAP(:,:,itime));
    elseif isfield(AQ,'SMOS')
      A=squeeze(AQ.SMOS(:,:,itime));
    end
    
    I=find(~isnan(A));
    smm(I)=smm(I)+A(I);
    ncc(I)=ncc(I)+1;
  end
end
I=find(ncc>0);
S=A*nan;
S(I)=smm(I)./ncc(I);
S(S==0)=nan;

c1=30;
c2=35;
nint = 200;
cmp = colormap(parula(nint));
cnt = (c1:(c2-c1)/nint:c2);
%CMP = colormap_sclr3(nint,c1,c2);
%cmp = CMP.colormap;
%cnt = CMP.intervals;
%nint=length(cmp);

% Greenland
xlim1 = 450;
xlim2 = 1100;
ylim1 = 300;
ylim2 = 1100;

hmsk=HH;
hmsk(HH<0)=nan;


figure(1); clf;
pcolor(I_hycom,J_hycom,S); shading flat;
hold on;
%contour(HH,[0 0],'k','Linewidth',1);
contour(HH,[-500 -500],'Color',[0.7 0.7 0.7],'Linewidth',0.5);
colormap(cmp);
caxis([c1 c2]);
axis('equal');
set(gca,'xlim',...
	[xlim1 xlim2],...
	'ylim',[ylim1 ylim2],...
	'xtick',[],'ytick',[],...
	'Color',[0.5 0.5 0.5]);
clr=[0.9 0.9 0.9];
plot_gridlines(45,10,1,clr,LON,LAT);


freezeColors;
pcolor(hmsk); shading flat;
colormap([0 0 0]);
freezeColors;

contour(I_hycom,J_hycom,S,[30:0.5:36],'k');
contour(I_hycom,J_hycom,S,[33 33],'k','Linewidth',1.6);

colormap(cmp);
hb = colorbar;
set(hb,'Position',[0.85 0.2 0.015 0.65],...
       'Fontsize',16,...
       'TickLength',0.02);

stl=sprintf('%s S, %i-%i, Mo:%i-%i',satl,iyr1,iyr2,mplt(1),mplt(end));
title(stl,'Fontsize',12,'Interpreter','none');

btx='plot_AquariusS.m';
bottom_text(btx,'pwd',1);

