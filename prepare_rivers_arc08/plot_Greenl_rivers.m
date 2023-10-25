% Monthly Greenland runoff prepared
% frmo Bumber's data: read_Greenland_rivers.m
% Note, Bamber's data are TOTAL Freshwater Flux
% that includes Solid Ice discharge, tundra runoff and
% Ice Sheet runoff
% Note: Greenalnd runoff after 2010 is extrapolated
% based on the trend estimate
% see: extrapolate_Greenl_rivers.m
%
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps
startup;

clear all
close

regn='ARCc0.08';
expt='050';  % hycom experiment *10
ntopo=9;      % HYCOM topo version
YearGr=2006;  % Gr. river runoff river
YR1 = 1958;
YR2 = 2015;
s_fig=logical(0);


PTH.data='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/force/rivers/';
PTH.topo='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
PTH.river='/Net/ocean/ddmitry/arctic_AOregimes/data/Greenland_rivers/';
PTH.fig = '/Net/tholia/ddmitry/hycom/ARCc0.08/051/fig_river/';
PTH.mat = '/Net/tholia/ddmitry/hycom/ARCc0.08/data_mat/';

if ~exist(PTH.fig,'dir')
  system(sprintf('mkdir -pv %s',PTH.fig));
end


mday=[31;28;31;30;31;30;31;31;30;31;30;31];
if mod(YearGr,4)==0,
  mday(2)=29;
end


ftopo = sprintf('%s/depth_%s_09.nc',PTH.topo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');


% Greenland:
fgrgr=sprintf('%sGreenland_grid.mat',PTH.river);
fgrrv=sprintf('%sGreenland_runoff_monthly-%4.4i.mat',PTH.river,YearGr);
% Greenland runoff and grid:
load(fgrrv);
load(fgrgr);
LN=GRgrd.LN;
LT=GRgrd.LT;
clear GRgrd
[mg,ng]=size(LN);

% COnvert LON/LAT to I,J
% 
fmtind=sprintf('%sGreen2HYCOMind.mat',PTH.mat);
get_ind=0;
if get_ind>0
  IG=[];
  JG=[];
  s = matlabpool('size');
  if s==0
    matlabpool open 12;
  end
  
  for ii=1:ng
    fprintf('ii=%i\n',ii);
    LNa=LN(:,ii);
    LTa=LT(:,ii);
    parfor jj=1:mg
%      ln0=LN(jj,ii);
%      lt0=LT(jj,ii);
      ln0=LNa(jj);
      lt0=LTa(jj);
      D=distance_spheric_coord(LAT,LON,lt0,ln0);
      [j0,i0]=find(D==min(min(D)));
      IGa(jj,1)=i0;
      JGa(jj,1)=j0;
    end
    IG(:,ii)=IGa;
    JG(:,ii)=JGa;
  end
  fprintf('Saving %s\n',fmtind);
  save(fmtind,'IG','JG');
  
    matlabpool close;
  
else
  load(fmtind);
end


%HG=interp2(LON,LAT,HH,LN,LT);

%fplt=logical(0);
% Convert km3/mo to m3/sec:
km3mo=1; % 1km3/mo
ndays=30.4167;
mSv=km3mo*1e9/(3600*24*ndays)*1e-3;  % mSv


c1=0;
c2=0.05;
nint=100;
cnt=(c1:(c2-c1)/nint:c2);  % construct array of intervals
Rr=-1;                        % colors of max intensity, clockwise rotation 
C0=[-1,0,1];                 % starting point, red dimension is Z axis
Cend=[1, 1, 0.6];
%cmp = colormap_spiral(nint,'C0',C0,'Rr',Rr,'Cend',Cend);
cmp = colormap_red(100);


nm=7;
A=GR(nm).runoff;

A(A==0)=nan;

[J,I]=find(~isnan(A));
In=find(~isnan(A));
clear C;
for jj=1:length(J);
  j0=J(jj);
  i0=I(jj);
  aa=A(j0,i0);
  ik=max(find(cnt<=aa));
  if ik>nint, ik=nint; end;
  if isempty(ik); ik=1; end;
  clr=cmp(ik,:);
  C(jj,:)=clr;
end
S=4;

hmm=HH;
hmm(hmm>0)=nan;
hmm(hmm<0)=-2;
hmm(1500:end,:)=nan;
hmm(:,1:400)=nan;
hmm(:,1000:end)=nan;

figure(1);
clf
pcolor(hmm); shading flat;
colormap([0 0 0]);
set(gca,'Color',[0.7 0.7 0.7]);
hold on;
freezeColors;
%scatter(I,J,16,C,'filled');
scatter(IG(In),JG(In),12,C,'filled');
%pcolor(A); shading flat;
%contour(HH,[0 0]
axis('equal');
hold on;
contour(LON,[-100:10:100],'Color',[0.8 0.8 0.8]);
contour(LAT,[40:10:90],'Color',[0.8 0.8 0.8]);
caxis([c1 c2]);
colormap(cmp);


%set(gca,'color',[1 1 1]);
%set(gca,'xlim',[0 300],'ylim',[0 560],'xtick',[],'ytick',[]);
set(gca,'xlim',[480 1000],'ylim',[380 1100],'xtick',[],'ytick',[]);
set(gcf,'Color',[1 1 1],'InvertHardCopy','off');
stt=sprintf('Total FW Flux, km3/mo, %4.4i %2.2i',YearGr,nm);
title(stt,'Fontsize',14);
%colorbar;
% Colorbar
hght=[];
lngth=[];
mint=10;
mbx=mint;
fsz=12;
bxc='k';
posc=[0.76 0.1 0.8 0.06];
aend=0;
[az,axc]  = colorbar_vert(cmp,cnt,hght,lngth,mint,fsz,bxc,posc,mbx,aend);
 
%keyboard

if s_fig
  fgout=sprintf('%sGrRunoff_%4.4i_%2.2i',PTH.fig,YearGr,nm);
  fprintf('Saving %s\n',fgout);
  print('-dpng','-r250',fgout);
end

%AOMIP river runoff: m3/s
Yen=[6038.57, 6022.93, 5983.87,6001.32, 27533.5, 77386.6, 26586.8, ...
      17485.4, 16896.1, 13969.2, 6855.90, 5839.77];
Len=[2783.04, 2136.78, 1651.78, 1350.20, 6235.90, 73917.4, 39683.4, ...
      27340.0, 24126.0, 13771.7, 3502.06, 2927.86];
Ob=[4986.72, 4120.23, 3635.81, 3697.96, 15122.7, 36715.9, 31694.8, ...
    23122.7, 14747.4, 11000.7, 6695.36, 5734.18]; %Ob+Pur
Mck=[3814.35, 3605.90, 3349.00, 3384.45, 13218.8, 21413.0, 17854.7, ...
     13984.2, 11268.7, 9038.15, 4756.00, 3595.00];
Pch=[959.400, 773.800, 695.400, 950.200, 15502.6, 17126.4, 5534.20, ...
     3227.80, 3917.00, 4197.60, 1894.40, 1277.40];

% m3/s -> km3/mo:
Ryen = Yen'*3600*24.*mday*1e-9;
Rlen = Len'*3600*24.*mday*1e-9;
Rob  = Ob'*3600*24.*mday*1e-9;
Rmck = Mck'*3600*24.*mday*1e-9;
Rpch = Pch'*3600*24.*mday*1e-9;
Rtot=Ryen+Rlen+Rob+Rmck+Rpch;


% =============
% Plot mo clim:
% =============
pclim=0;
if pclim>0
  YY=[1990;1998;2007];
  ny=length(YY);
  clear Rr
  for iy=1:ny
    Yrf=YY(iy);
    fgrrv=sprintf('%sGreenland_runoff_monthly-%4.4i.mat',PTH.river,Yrf);
    fprintf('%s\n',fgrrv);
    load(fgrrv);  
    for ii=1:12
      A=GR(ii).runoff;
      Rr(ii,iy)=nansum(nansum(A)); % km3/mo
    end
    LGD{iy}=sprintf('%4.4i',Yrf);
  end;
  CLR=[0,0,1; 0,1,0; 1,0,0; 0,0,0];

  nf=2;
  figure(nf); clf;
  hold on
  for iy=1:ny
    clr=CLR(iy,:);
    plot(Rr(:,iy),'Color',clr,'linewidth',2);
  end;
  %plot(Rtot,'k','linewidth',2);
  %LGD{end+1}='RArct';
  legend(LGD,'Fontsize',12);

  set(gca,'tickdir','out','xlim',[1 12],'ylim',[0 300],'Fontsize',12);
  set(gca,'xtick',[1:12],'ytick',[0:25:300],'xgrid','on','ygrid','on');
  set(gcf,'Color',[1 1 1],'InvertHardCopy','off');
  xlabel('Months');
  ylabel('Total FW Flux, km^3/mo');
  title('FW Flux, Greenland');


  if s_fig
    fgout=sprintf('%smonthlyGrRunoff',PTH.fig);
    fprintf('Saving %s\n',fgout);
    print('-dpng','-r150',fgout);
  end
end




% Plot interannual variability:
clear Rt Rym
cc=0;
for iy=YR1:YR2
  Yrf=iy;
  fgrrv=sprintf('%sGreenland_runoff_monthly-%4.4i.mat',PTH.river,Yrf);
  fprintf('%s\n',fgrrv);
  load(fgrrv);  
  A=GR(7).runoff;
  cc=cc+1;
  Rt(cc,1)=nansum(nansum(A)); % km3/mo
  
  for im=1:12
    A=GR(im).runoff;
    Rym(cc,im)=nansum(nansum(A));
  end
  
end;

tm=[YR1:YR2];
mo=[1:12];


Rann = sum(Rym,2);

figure(4); clf;
pcolor(mo,tm,Rym); shading flat;
colorbar
title('Total Greenland runoff, km^3/mo');


figure(3); clf;
axes('Position',[0.08 0.55 0.87 0.4]);
plot(tm,Rann,'Linewidth',2);
set(gca,'tickdir','out',...
	'xlim',[tm(1) tm(end)],...
	'ylim',[699 1253],...
	'fontsize',12,...
        'xtick',[1950:5:2010],...
	'ytick',[600:50:3000],...
	'xgrid','on',...
	'ygrid','on');
title('Annual Greenland Runoff, km^3/mo');

txtb = 'plot_Greenland_rivers.m';
bottom_text(txtb,'pwd',1,'position',[0.08 0.4 0.8 0.1]);

%plot(tm,Rt,'Linewidth',2);
%set(gca,'tickdir','out','xlim',[tm(1) tm(end)],'ylim',[140 280],'fontsize',12);
%set(gca,'xtick',[150:5:2010],'ytick',[100:25:300],'xgrid','on','ygrid','on');
%set(gcf,'Color',[1 1 1],'InvertHardCopy','off');
%title('July Greenland FW Flux, km^3/mo');

if s_fig
%  fgout=sprintf('%sJulyGrRunoff',PTH.fig);
  fgout=sprintf('%sAnnualGrRunoff',PTH.fig);
  fprintf('Saving %s\n',fgout);
  print('-dpng','-r150',fgout);
end

