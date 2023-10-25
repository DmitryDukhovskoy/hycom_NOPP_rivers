% Plot T/S variability in Irminger Sea near-surface
% http://hadobs.metoffice.com/en4/index.html
%The EN4 dataset consists of two products:
%Observed subsurface ocean temperature and 
% salinity profiles with data quality information, and,
%Objective analyses formed from the profile data with uncertainty estimates.
%Data are available from 1900 to the present and 
% there are separate files for each month.
% Please read 'Good, S. A., M. J. Martin and N. A. Rayner, 2013. EN4: 
% quality controlled ocean temperature and salinity profiles and 
% monthly objective analyses with uncertainty estimates, 
% Journal of Geophysical Research: Oceans, 118, 6704-6716, 
% doi:10.1002/2013JC009067' for details of how the dataset was constructed.
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

s_fig = 1; 
s_mat = 1;
s_calc = 0; % derive depths

%pthdat = '/Net/yucatan/tachanat/ocean_analysis/EN4/EN4_extract/';
pthdat = '/Net/kronos/ddmitry/EN4/';
pthmat = '/Net/ocean/ddmitry/vector_winds/dataM/'; 
pthmat = '/nexsan/people/ddmitry/data_mat/';
pthfig = '/Net/ocean/ddmitry/vector_winds/fig_EN4/';
en4v   = 'EN.4.2.0.f.analysis.g10'; % EN4 version
txtb   = 'ts_IrmingerSea.m';
fmat   = sprintf('%sts_Irminger_EN4.mat',pthmat);

%YC1=1993;
%YC2=1993;
yr1=2000;
yr2=2015;
%fclim = sprintf('%sEN4_climatology%i_%i.mat',pthmat,YC1,YC2);
%fanom = sprintf('%sEN4_climatology%i_%i.mat',pthmat,yr1,yr2);

fnm = sprintf('%s%s.201302.nc',pthdat,en4v);

S = double(nc_varget(fnm,'salinity'));
S = squeeze(S(1,1,:,:));

% Find sections:
LON = nc_varget(fnm,'lon');
LAT = nc_varget(fnm,'lat');
ZZ  = -1*nc_varget(fnm,'depth');
I = find(LON>180);
LON(I)=LON(I)-360;

% Depth in the central Irminger Sea: 
i1=max(find(LON<=-43));
i2=max(find(LON<=-37));
j1=min(find(LAT>=58));
j2=max(find(LAT<=61));

SGMD.Name   = sprintf('IrmingerSea, T/S spatially average');
SGMD.Code   = 'anls_mtlb_utils/hycom_NOPP/EN4/ts_IrmingerSea.m';
SGMD.IJ     = [i1,j1;i2,j2];
SGMD.YCoord = [LAT(j1),LAT(j2)];
SGMD.XCoord = [LON(i1),LON(i2)];


cc=0;
for ii=yr1:yr2
  for im=1:12
    if im>2 & im<11, continue; end
    cc=cc+1;
    YRPLT(cc,1)=ii;
    YRPLT(cc,2)=im;
  end
end

if s_calc>0
  nrc=cc;
  cc=0;
  for ik=1:nrc
    YR=YRPLT(ik,1);
    IM=YRPLT(ik,2);
    fnm = sprintf('%s%s.%4.4i%2.2i.nc',pthdat,en4v,YR,IM);

    fprintf('EN4: Reading %i/%i, %s\n',YR,IM,fnm);

    S = double(nc_varget(fnm,'salinity'));
    T = double(nc_varget(fnm,'temperature'))-273.15; % K -> C
    cc=cc+1;

    S = squeeze(S(1,:,:,:));
    T = squeeze(T(1,:,:,:));
    [ll,mm,nn]=size(S);

    clear SS TT
    if i2<i1,
      S1=S(:,j1:j2,i1:end);
      S2=S(:,j1:j2,1:i2);
      n1=size(S1,3);
      n2=size(S2,3);
      SS=S1;
      for in=1:n2
	n1=n1+1;
	SS(:,:,n1)=S2(:,:,in);
      end
      T1=T(:,j1:j2,i1:end);
      T2=T(:,j1:j2,1:i2);
      n1=size(T1,3);
      n2=size(T2,3);
      TT=T1;
      for in=1:n2
	n1=n1+1;
	TT(:,:,n1)=T2(:,:,in);
      end
    else
      TT = T(:,j1:j2,i1:i2);
      SS = S(:,j1:j2,i1:i2);
    end

    dmm=squeeze(SS(1,:,:));
    I=find(~isnan(dmm));
    Sav = squeeze(nanmean(SS(:,I),2));
    Tav = squeeze(nanmean(TT(:,I),2));
  
    nl = size(Sav,1);
    SGMD.TM(cc,1) = datenum(YR,IM,1);
    SGMD.Tav(cc,1:nl) = Tav;
    SGMD.Sav(cc,1:nl) = Sav;
    SGMD.ZZ   = ZZ;
    fprintf('Surf Tav = %6.2f, Sav = %6.2f\n',Tav(1),Sav(1));

  end % Time

  if s_mat>0
    fprintf('Saving %s\n',fmat);
    save(fmat,'SGMD');
  end
else
  fprintf('Loading %s\n',fmat);
  load(fmat);
end

TM = SGMD.TM;
DV = datevec(TM);
Tav = SGMD.Tav;
Sav = SGMD.Sav;
ZZ  = SGMD.ZZ;
DV = datevec(TM);
YR1 = DV(1,1);
yrs = [1:length(TM)]./4+YR1;

iz1 = 1;
iz2 = 20;
iz3 = 28;

s1 = Sav(:,iz1);
s2 = Sav(:,iz2);
s3 = Sav(:,iz3);
z1 = abs(round(ZZ(iz1)));
z2 = abs(round(ZZ(iz2)));
z3 = abs(round(ZZ(iz3)));

nl = size(Sav,2);
for ik=1:nl
  s10(ik,1)=prctile(Sav(:,ik),10);
  s90(ik,1)=prctile(Sav(:,ik),90);
end
sM=mean(Sav,1)';


figure(1); clf;
axes('position',[0.08 0.08 0.25 0.8]);
plot(sM,ZZ,'r-');
hold on
plot(s10,ZZ,'b-');
plot(s90,ZZ,'b-');
set(gca,'xlim',[34.63 35.01],...
	'ylim',[-1330 0]);

title('EN4, S, 2000-2015');

axes('position',[0.39 0.5 0.58 0.35]);
plot(yrs,s1,'r-');
hold on;
plot(yrs,s2,'b-');
plot(yrs,s3,'g-');
set(gca,'xlim',[2000 2016]);

ttl=sprintf('EN4, S %im (red), %im (blue), %i(gr)',z1,z2,z3); 
title(ttl);

axes('position',[0.4 0.08 0.25 0.2]);
hist(s1,30);
title(sprintf('S, %i m',z1));
set(gca,'xlim',[34.6 35.2]);

axes('position',[0.73 0.08 0.25 0.2]);
hist(s2,30);
title(sprintf('S, %i m',z2));
set(gca,'xlim',[34.6 35.2]);

%bottom_text(txtb,'pwd');

if s_fig>0
  fgnm=sprintf('%sS_EN4_IrmingerSea',pthfig);
  fprintf('Saving %s\n',fgnm);
  print('-dpng','-r0',fgnm);
end

% T 
t1 = Tav(:,iz1);
t2 = Tav(:,iz2);
t3 = Tav(:,iz3);

nl = size(Tav,2);
for ik=1:nl
  t10(ik,1)=prctile(Tav(:,ik),10);
  t90(ik,1)=prctile(Tav(:,ik),90);
end
tM=mean(Tav,1)';

figure(2); clf;
axes('position',[0.08 0.08 0.25 0.8]);
plot(tM,ZZ,'r-');
hold on
plot(t10,ZZ,'b-');
plot(t90,ZZ,'b-');
set(gca,'xlim',[0.93 9],...
	'ylim',[-1330 0]);

ttl=sprintf('EN4, T %im (red), %im (blue), %i(gr)',z1,z2,z3); 
title(sprintf('T, %i m',z1));

axes('position',[0.39 0.5 0.58 0.35]);
plot(yrs,t1,'r-');
hold on;
plot(yrs,t2,'b-');
plot(yrs,t3,'g-');
set(gca,'xlim',[2000 2016]);

title(ttl); 

axes('position',[0.4 0.08 0.25 0.2]);
hist(t1,30);
title(sprintf('T, %i m',z1));

axes('position',[0.73 0.08 0.25 0.2]);
hist(t2,30);
title(sprintf('T, %i m',z2));


%bottom_text(txtb,'pwd');

if s_fig>0
  fgnm=sprintf('%sT_EN4_IrmingerSea',pthfig);
  fprintf('Saving %s\n',fgnm);
  print('-dpng','-f2','-r0',fgnm);
end

