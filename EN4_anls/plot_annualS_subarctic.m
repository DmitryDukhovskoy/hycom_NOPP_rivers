% annual mean S in North Atlantic
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

%pthdat = '/Net/yucatan/tachanat/ocean_analysis/EN4/EN4_extract/';
pthdat = '/Net/kronos/ddmitry/EN4/';
%pthmat = '/Net/ocean/ddmitry/vector_winds/dataM/'; 
pthmat = '/nexsan/people/ddmitry/data_mat/';
pthmat2= '/Net/mars/ddmitry/hycom/GLBb0.08/data_mat/';
pthfig = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/fig_EN4/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
en4v   = 'EN.4.2.0.f.analysis.g10'; % EN4 version
txtb   = 'annualS_subarctic.m';

LR=[0, -51;...
    -51,-151;...
    -151,-301];
nlrs=size(LR,1);

%ftopo = sprintf('%sGLBb_T07_subset_Natl.mat',pthmat2);
%fprintf('Loading topo %s\n',ftopo);
%TH = load(ftopo);
%HH = TH.HH;
%LONH = TH.LON;
%LATH = TH.LAT;

ftopo = sprintf('%s/depth_ARCc0.08_11.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
LONH = nc_varget(ftopo,'Longitude');
LATH = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LONH);


%Sref = 35.5;
%Sref = 35; 
yr1=1983;
yr2=2016;
fmat = sprintf('%sannualS_subarctic_EN4.mat',pthmat);

fprintf('Loading %s\n',fmat);
load(fmat);

% average S by regions
BX = sub_define_boxes(HH,LONH,LATH,0);

[XX,YY] = meshgrid((1:nn),(1:mm));
for ib=1:5
  iBG = BX(ib).IJ;
  INP = inpolygon(XX,YY,iBG(:,1),iBG(:,2));
  IN = find(INP==1);
  BX(ib).IN = IN;
end


xmn = min(SANN(1).LON);
if xmn>=0;
  lnn=SANN(1).LON;
  lnn(lnn>180)=lnn(lnn>180)-360;
  SANN(1).LON=lnn;
end

LON=SANN(1).LON;
LAT=SANN(1).LAT;
[LN,LT] = meshgrid(SANN(1).LON,SANN(2).LAT);
[DX,DY] = sub_dx_dy(LN,LT);
Acell = DX.*DY*1e-6; % km2

for ib=1:5
  Xb=BX(ib).XY(:,1);
  Yb=BX(ib).XY(:,2);
  xmn=min(Xb);
  xmx=max(Xb);
  
  if sign(xmn)==sign(xmx)
    inp = inpolygon(LN,LT,Xb,Yb);
    IN  = find(inp==1);
  else
    nL=length(LON);
    nT=length(LAT);
    x1=min(Xb);
    x2=max(Xb);
    y1=min(Yb);
    y2=max(Yb);
    i1=max(find(LON<=x1));
    i2=max(find(LON<=x2 & LON>0));
    j1=max(find(LAT<=y1));
    j2=max(find(LAT<=y2));
    I1=[i1:nL];
    I2=[1:i2];
    J1=[j1:j2];
    II=[I1,I2];
    cx=0;
    for ii=II
      for jj=J1
	cx=cx+1;
	I=sub2ind(size(LN),jj,ii);
	IN(cx,1)=I;
      end
    end
  end
  BX(ib).IN=IN;
end

f_chck=0;
if f_chck==1
  x=SANN(1).LON;
  y=SANN(1).LAT;
  i=find(x<=0);
  x(i)=x(i)+360;
%  x(end)=360;
  S=squeeze(SANN(1).annualS(ilv,:,:));
  pcolor(x,y,S); shading flat;
  hold on
  
  I=find(LN<=0);
  LN(I)=LN(I)+360;
  for ib=1:5
    IN=BX(ib).IN;
    plot(LN(IN),LT(IN),'k.');
  end
end

nyrs=length(SANN);
for ir=1:nyrs  
  for ilv=1:3
    S=squeeze(SANN(ir).annualS(ilv,:,:));
    
    for ib=1:5
      IN=BX(ib).IN;
      sav=nanmean(S(IN));
      SAV(ib).S(ilv,ir)=sav;
    end
    
  end
end

POS(1,:)=[0.08 0.67 0.85 0.25];
POS(2,:)=[0.08 0.35 0.85 0.25];
POS(3,:)=[0.08 0.08 0.85 0.25];
btx = 'plot_annualS_subarctic.m';
YRS=[yr1:yr2];
ryr=1993;
iyr=find(YRS==ryr);
yrs=[ryr:2016];
for ib=1:5
  nm=BX(ib).Name;
  figure(ib); clf
  axes('Position',[0.08 0.4 0.85 0.5]);
  ss=SAV(ib).S;
  clear dS
  for ilv=1:3
    s90=mean(ss(ilv,1:iyr-1)); % mean 1983-1992
    ds=ss(ilv,:)-s90;
    dS(:,ilv)=ds(iyr:end)';
  end
  hb=bar(yrs,dS);
  stl=sprintf('EN4-4.2.1, %s, dltS wrt S(1983-1992)=%4.2f',nm,s90);
  title(stl);
  legend('50m','150m','300m');
  
  set(gca,'tickdir','out',...
	  'xlim',[1993-0.6 2016+0.6],...
	  'xtick',[1990:2020]);
  
  bottom_text(btx,'pwd',1,'position',[0.05 0.3 0.6 0.06]);
  
end

% Estimate FW content change
% using mean S before 1993 as a reference
% Baffin Bay only, as it only has -dltS
ib=5;
nm=BX(ib).Name;
ss=SAV(ib).S;
clear dS
IN=BX(ib).IN;
Areg = sum(Acell(IN))*1e6; % m2
for ilv=1:3
  s90=mean(ss(ilv,1:iyr-1));
  ds=ss(ilv,iyr:end)-s90;
  dz = abs(LR(ilv,2)-LR(ilv,1));
  v2 = Areg*dz;  %m3
  vfw = v2-v2*(s90+ds)/s90;
end

figure(10); clf
axes('Position',[0.08 0.4 0.85 0.5]);
hb=bar(vfw./Areg); 

  
