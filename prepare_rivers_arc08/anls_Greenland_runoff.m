% Plot Bamber river data
% km3/mo
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
startup

close all
clear

s_fig=0;

hg=2^100; 

PTH.data='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/force/rivers/';
PTH.topo='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
PTH.river='/Net/ocean/ddmitry/arctic_AOregimes/data/Greenland_rivers/';
pthmat  = '/Net/tholia/ddmitry/hycom/ARCc0.08/data_mat/';
pthfig  = '/Net/tholia/ddmitry/hycom/ARCc0.08/fig_Greenland_river/';

YearGr=2004;
mday=[31;28;31;30;31;30;31;31;30;31;30;31];
if mod(YearGr,4)==0,
  mday(2)=29;
end
fgrgr=sprintf('%sGreenland_grid.mat',PTH.river);
load(fgrgr);
LN=GRgrd.LN;
LT=GRgrd.LT;

% Get rivers by regions:
% Approximate regions Baffin Bay, Labr. Sea, Nordic Seas:
% [I, J] indices on Bamber;s grid
rBB=[1,200; 100,200; 100,450; 1,450];
rLS=[20,200; 100,200; 100,2; 20,2];
rNS=[240,500; 300,500; 300,220; 120,220; 120,450; 238,501]; % Nord. Seas
rAO=[1,450; 120,450; 238,501; 240,500; 240,560; 1,560];   % Arctic ocean
rIS=[100,2; 100,220; 300,220; 300,2];  % Irminger Sea

rBB(end+1,:)=rBB(1,:);
rLS(end+1,:)=rLS(1,:);
rNS(end+1,:)=rNS(1,:);
rAO(end+1,:)=rAO(1,:);
rIS(end+1,:)=rIS(1,:);

[a1,a2]=size(LN);
[X,Y]=meshgrid([1:a2],[1:a1]);
in=inpolygon(X,Y,rBB(:,1),rBB(:,2)); % Baffin
iBB=find(in==1);
in=inpolygon(X,Y,rLS(:,1),rLS(:,2)); % Labrador
iLS=find(in==1);
in=inpolygon(X,Y,rNS(:,1),rNS(:,2)); % Nordic
iNS=find(in==1);
in=inpolygon(X,Y,rAO(:,1),rAO(:,2)); % AO
iAO=find(in==1);
in=inpolygon(X,Y,rIS(:,1),rIS(:,2)); % Irminger
iIS=find(in==1);
IR(1).Name='BB';
IR(1).I=iBB;
IR(2).Name='LS';
IR(2).I=iLS;
IR(3).Name='NS';
IR(3).I=iNS;
IR(4).Name='AO';
IR(4).I=iAO;
IR(5).Name='IS';
IR(5).I=iIS;

% Greenland:
rv61_90=876; % mean total runoff over 1961-1990
             % estimate as mean(rv), for YR1=1961, YR2=1990;
nyr=0;
YR1=1958;
%YR2=1989;
%YR1=1990;
YR2=2010;
clear rv Rr
for iyear=YR1:YR2
  YearGr=iyear;
  fgrgr=sprintf('%sGreenland_grid.mat',PTH.river);
  fgrrv=sprintf('%sGreenland_runoff_monthly-%4.4i.mat',PTH.river,YearGr);
% Greenland runoff and grid:
  load(fgrrv); % -> GR struct with monthly river maps
  load(fgrgr);
  Xgr=GRgrd.LN;
  Ygr=GRgrd.LT;
  clear GRgrd

  nyr=nyr+1;
% Get Greenland total monthy runoff:
  for im=1:12
    a=GR(im).runoff; % km3/mo
    rv(im,nyr)=nansum(nansum(a));
    if im==7 & iyear==2004,
      a7=a;
    end
%
% Regions:
    for jr=1:5
      I=IR(jr).I;
      Rr(jr).Region=IR(jr).Name;
      Rr(jr).river(im,nyr)=nansum(a(I));
    end
    
  end
end;


btmtxt='anls_Greenland_runoff.m';

figure(10); clf; 
hold on;
pcolor(a7); shading flat;
axis('equal');
set(gca,'xlim',[-1 400],'ylim',[1 550]);
plot(rBB(:,1),rBB(:,2),'r.-','Linewidth',2);
plot(rLS(:,1),rLS(:,2),'g.-','Linewidth',2);
plot(rNS(:,1),rNS(:,2),'b.-','Linewidth',2);
plot(rAO(:,1),rAO(:,2),'k.-','Linewidth',2);
plot(rIS(:,1),rIS(:,2),'m.-','Linewidth',2);
%set(H,'Position',[0.2 0.5]);
bottom_text(btmtxt,'pwd',1);

if s_fig>0
  fgnm=sprintf('%sanls_GrRunoff_regions',pthfig);
  fprintf('Saving %s\n',fgnm);
  print('-dpng','-r300',fgnm);
end



atot=sum(rv,1); % annual runoff km3/yr
Rmin=min(rv,[],2); % min runoff over N years
Rmax=max(rv,[],2);

% Analysis of cumulative river runoff
% using all data: need to change YR1 to 1958 above
yr=[YR1:YR2];
%a1=atot(1:32);
%a2=atot(33:end);
%cs=sum(a2-mean(a1));


% Add new Bamber's data 
% for comparison:
f_newdata=1;
if f_newdata==1
  fmat = sprintf('%sGreenlFWFv3.mat',pthmat);
  fprintf('Loading %s\n',fmat);
  load(fmat)
  
  YRn = GrF.Years;
  FWFn = GrF.FWFtotal_km3yr;
  FWF  = GrF.TotalFWF_km3mo;
else
  YRn =yr;
  FWFn = [];
end

% Plot July FWF from old and new data sets
nrc = length(FWF);
nyr = nrc/12;
FTMn = reshape(FWF,[12, nyr]);
imo = 3;
fnew = FTMn(imo,:);
fold = rv(imo,:);

yr1 = min([YRn(1), yr(1)]);
yr2 = max([YRn(end), yr(end)]);

figure(5); clf;
axes('Position',[0.08 0.5 0.85 0.43]);
plot(yr,fold);
hold on
plot(YRn,fnew,'r');
set(gca,'tickdir','out',...
	'xlim',[yr1 yr2],...
	'xtick',[1950:5:2016],...
	'xgrid','on',...
	'ygrid','on');
stl=sprintf('M=%i, GreenData Old (b) & New(r), Total FWF, km3/mo',imo);
title(stl);

btx = 'anls_Greenland_runoff.m';
bottom_text(btx,'pwd',1,'Position',[0.03 0.4 0.9 0.1]);


  
  
figure(1); clf;
axes('Position',[0.09 0.12 0.85 0.35]);
plot(yr,atot);
hold on
plot(YRn,FWFn,'r');
set(gca,'tickdir','out',...
	'xlim',[yr1 yr2],...
	'xtick',[1950:5:2016],...
	'xgrid','on',...
	'ygrid','on');
title('Data Old (b) and New(r), Total Greenl. FW Flux, km3/yr');

btx = 'anls_Greenland_runoff.m';
bottom_text(btx,'pwd',1);

%
%  Rivers by regions:
figure(2); clf;
hold on;

CLR=[0,0,0.6; 0,0.8,1; 0,1,0.4; 1,0.8,0; 1,0,0; 0.8,0.3,0];
yr=[YR1:YR2];
nyr=length(yr);
for jr=1:5  % regions
  rr=Rr(jr).river;
  rann=sum(rr,1);  % annual runoff in a region
  GrRg(jr,1:nyr)=rann;  % group by regions
  Rtot(jr)=sum(rann)*1e9;  % m3 - total runoff over all years, regions
  clr=CLR(jr,:);
  plot(yr,rann,'Color',clr,'Linewidth',2.5);
  ylabel('km3/yr');
  xlabel('years');
end
clr=CLR(6,:);
plot(yr,atot,'k-','Color',clr,'Linewidth',2);
legend('Baff','Labr','Nord','AO','Irm','total','Location','BestOutside');
title('Annual river runoff');
bottom_text(btmtxt,'pwd',1);
if s_fig>0
  fgnm=sprintf('%sanls_GrRunoff_annual_regions',pthfig);
  fprintf('Saving %s\n',fgnm);
  print('-dpng','-r300',fgnm);
end

% =======================
% Bar diagram of annual runoff by regions
% =======================
cmp=CLR(1:5,:);
figure(3); clf;
axes('Position',[0.08 0.35 0.9 0.58]);
hb=bar(yr,GrRg','stacked');
colormap(cmp);
hold on;
set(gca,'tickdir','out',...
	'xlim',[1989 2011],...
	'xtick',[1990:2010],...
	'ylim',[0 1300],...
	'ytick',[0:200:1300],...
	'yminortick','on',...
	'fontsize',12);
hg=legend('Baff','Labr','Nord','AO','Irm');
set(hg,'Position',[0.5 0.05 0.1 0.2],'fontsize',10);
plot(yr,atot,'k-','Color',[0 0 0],'linewidth',2);
title('Annual Greenland Freshwater flux (km3) in 5 regions & total');
bottom_text(btmtxt,'pwd',1);
if s_fig>0
  fgnm=sprintf('%sanls_GreenlRunoff_barAnnual',pthfig);
  fprintf('Saving %s\n',fgnm);
  print('-dpng','-r300',fgnm);
  print('-depsc2',fgnm);
end





% Runoff over all years;
Ftot=sum(atot)*1e9;  % m3
fprintf('Total runoff over %i years = %9.3fe12 m3\n',...
	nyr, Ftot*1e-12);
fprintf('  Total Baffin runoff over %iyrs = %9.3fe12 m3: \n',...
	nyr,Rtot(1)*1e-12);
fprintf('  Total Labrador runoff over %iyrs = %9.3fe12 m3: \n',...
	nyr,Rtot(2)*1e-12);
fprintf('  Total Nordic runoff over %iyrs = %9.3fe12 m3: \n',...
	nyr,Rtot(3)*1e-12);
fprintf('  Total AO runoff over %iyrs = %9.3fe12 m3: \n',...
	nyr,Rtot(4)*1e-12);
fprintf('  Total Irm runoff over %iyrs = %9.3fe12 m3: \n',...
	nyr,Rtot(5)*1e-12);





% Plot monthly histograms 
Nyr=7;
Rtot=zeros(12,1);
figure(4); clf;
for jr=1:6
  if jr<6
    rname=IR(jr).Name;
    rr=Rr(jr).river;
    rm=mean(rr,2);
    Rtot=Rtot+rm;
  else
    rname='tot';
    rm=Rtot;
  end
  
  dy=0.24;
  dx=0.3;
  switch(jr)
   case(1)
    pp=[0.1 0.7 dx dy];
   case(2)
    pp=[0.5 0.7 dx dy];
   case(3)
    pp=[0.1 0.38 dx dy];
   case(4)
    pp=[0.5 0.38 dx dy];
   case(5)
    pp=[0.1 0.06 dx dy];
   case(6)
    pp=[0.5 0.06 dx dy];
  end
  
  axes('Position',pp);
  hb=bar(rm,0.9);
  hold on;
  set(hb,'FaceColor',[0.7 0.7 0.7]);
  if jr<6
    for ii=1:12
      rmin=min(rr(ii,:));
      rmax=max(rr(ii,:));
      plot([ii ii],[rmin rmax],'k-','Linewidth',3);
      yuplt=1.1*max(max(rr));
    end;
  else
    for ii=1:12
      rmin=Rmin(ii);
      rmax=Rmax(ii);
      plot([ii ii],[rmin rmax],'k-','Linewidth',3);
      yuplt=1.1*max(Rmax);
    end;
  end    
  ftot=sum(rm);
  ftxt=sprintf('%6.1f km3/yr',ftot);
  ttl=sprintf('%s, km3/mo',rname);
  title(ttl);
  set(gca,'xlim',[0 12.5],...
	  'ylim',[0 yuplt],...
	  'tickdir','out',...
	  'Fontsize',11);

  text(0.5,0.85*yuplt,ftxt);
  
end

bottom_text(btmtxt,'pwd',1);
if s_fig>0
  fgnm=sprintf('%sanls_GrRunoff_barSeasonal',pthfig);
  fprintf('Saving %s\n',fgnm);
  print('-depsc2',fgnm);
  print('-dpng','-r300',fgnm);
end

