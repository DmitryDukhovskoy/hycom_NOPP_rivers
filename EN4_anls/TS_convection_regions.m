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

s_fig = 0; 
s_mat = 0; % =0 - load existing mat file

%pthdat = '/Net/yucatan/tachanat/ocean_analysis/EN4/EN4_extract/';
pthdat = '/Net/kronos/ddmitry/EN4/';
%pthmat = '/Net/ocean/ddmitry/vector_winds/dataM/'; 
pthmat = '/nexsan/people/ddmitry/data_mat/';
pthfig = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/fig_EN4/';
en4v   = 'EN.4.2.0.f.analysis.g10'; % EN4 version
txtb   = 'ts_IrmingerSea.m';
fmat   = sprintf('%sTS_ConvectionSites_EN4.mat',pthmat);

%YC1=1993;
%YC2=1993;
yr1=1990;
yr2=2016;
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

TS = sub_regions;
nR = length(TS);


% Find indices:
[LN,LT] = meshgrid(LON,LAT);

for i=1:nR
  XY = TS(i).XY;
%  if i ~= 4
  pl = inpolygon(LN,LT,XY(:,1),XY(:,2));
  IN = find(pl==1);
  TS(i).IN = IN;
%  else

  TS(i).ZZ   = ZZ;
  TS(i).Code = 'anls_mtlb_utils/hycom_NOPP/EN4_anls/TS_convection_regions.m';
end


cc=0;
for ii=yr1:yr2
  for im=1:12
%    if im>2 & im<11, continue; end
    cc=cc+1;
    YRPLT(cc,1)=ii;
    YRPLT(cc,2)=im;
  end
end

if s_mat>0
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
    
    for i=1:nR
      IN = TS(i).IN;
      SS = S(:,IN);
      TT = T(:,IN);

      Sav = squeeze(nanmean(SS,2));
      Tav = squeeze(nanmean(TT,2));
      fprintf('Reg. %i: Surf Tav = %6.2f, Sav = %6.2f\n',i,Tav(1),Sav(1));
  
      nl = size(Sav,1);
      TS(i).TM(cc,1) = datenum(YR,IM,1);
      TS(i).Tav(cc,1:nl) = Tav;
      TS(i).Sav(cc,1:nl) = Sav;
    end
    

  end % Time

  if s_mat==1
    fprintf('Saving %s\n',fmat);
    save(fmat,'TS');
  end
  
else
  fprintf('Loading %s\n',fmat);
  load(fmat);
end


TM = TS(1).TM;
DV = datevec(TM);
YR1 = DV(1,1);
YR2 = DV(end,1);
yrs = [0:length(TM)-1]./12+YR1;
ZZ  = TS(1).ZZ;
lz = length(ZZ);
dz = abs(diff(ZZ));
dz(lz) = dz(lz-1);
clear ZZi
ZZi(1,1) = 0;
for kk=1:lz
  ZZi(kk+1) = -(abs(ZZ(kk))+abs(0.5*dz(kk)));
end
ZZi=ZZi(:);
Zdp=abs(diff(ZZi)); % layer thicknesses

iz1 = 1;
iz2 = 6; % upper 55 m
iz3 = 11; % upper 110 m
z1 = abs(round(ZZ(iz1)));
z2 = abs(round(ZZ(iz2)));
z3 = abs(round(ZZ(iz3)));

% For all regions, 
% find depth-averaged S over depths z1:z2
clear S1
for i=1:nR
  Tav = TS(i).Tav;
  Sav = TS(i).Sav;
  dmm = Sav(:,1:iz2);
  zp  = Zdp(1:iz2); % layer thickness
  [a1,a2] = size(dmm);
  for k=1:a1,
    spp = dmm(k,:);
    spp=spp(:);
    smm = spp'*zp/sum(zp);
    S1(i,k) = smm;
  end
end 

% Calculate yearly mean S:
[a1,a2]=size(S1);
Nyr=a2/12;
clear Sm
for i=1:nR
  dmm=S1(i,:);
  b=reshape(dmm,[12,Nyr]);
  Sm(i,:) = nanmean(b,1);
end



xlbl=[];
cll=0;
for yr=yrs(1):yrs(end);
  cll=cll+1;
  if mod(yr,5)==0
    xlbl{cll}=sprintf('%i',yr);
  else
    xlbl{cll}=' ';
  end
end


for ik=1:nR
  figure(ik); clf;
  axes('position',[0.1 0.5 0.85 0.4]);
  ss1=S1(ik,:);
  plot(yrs,ss1,'-');
  hold on;
  for jr=1:Nyr
    x1=YR1+jr-1;
    x2=x1+0.999;
    plot([x1 x2],[Sm(ik,jr) Sm(ik,jr)],'r-');
  end
  
  set(gca,'xlim',[min(yrs)-0.1 max(yrs)+0.1],...
	  'ylim',[0.999*min(ss1) 1.001*max(ss1)],...
	  'tickdir','out',...
	  'xtick',[yr1:yr2],...
	  'xticklabel',xlbl,...
          'xminortick','on',...
	  'xgrid','on',...
	  'ygrid','on');

  nm=TS(ik).Name;
  stl = sprintf('EN4, S avrg 0-%im, %i-%i, %s',...
		z2,YR1,YR2,nm);
  
  title(stl);
  btx = 'TS_convection_regions.m';
  bottom_text(btx,'pwd',1,'position',[0.08 0.3 0.7 0.2]);

  if s_fig>0
    fgnm=sprintf('%sS_EN4_%s_%im',pthfig,nm,z2);
    fprintf('Saving %s\n',fgnm);
    print('-dpng','-r250',fgnm);
  end

end



