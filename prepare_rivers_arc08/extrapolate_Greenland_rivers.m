% Extrapolate Greenland monthly runoff
% for the time period not covered in 
% bamber's data
%
% Bamber's data cover period until 2010
% after that, I'm adding a trend estimated
% from Bamber's total Greenland Mass Change
%  ~-290 GT/yr for 2010-2014
%

addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
startup;

clear all
close

f_save=1;  % save data in mat files

PTH.river='/Net/ocean/ddmitry/arctic_AOregimes/data/Greenland_rivers/';

YS=2011;
YE=2016;

% First esimate trend in annual runoff for 2003-2009
% to compare with trend in total mass change for 2003-2009
% which is ~235 Gt/year
% and derive the trend for 2010-2014
% total mass change is ~ -294 Gt/yr

YR1 = 1992;
YR2 = 2010;
nyr = 0;
clear rv Rr
for iyear=YR1:YR2
  YearGr=iyear;
  fgrgr=sprintf('%sGreenland_grid.mat',PTH.river);
  fgrrv=sprintf('%sGreenland_runoff_monthly-%4.4i.mat',PTH.river,YearGr);
% Greenland runoff and grid:
  fprintf('Loading %s\n',fgrrv);
  load(fgrrv); % -> GR struct with monthly river maps

  nyr=nyr+1;
% Get Greenland total monthy runoff:
  for im=1:12
    a=GR(im).runoff; % km3/mo
    rv(im,nyr)=nansum(nansum(a));
  end
end;

% Save last year:
GR0=GR;

% Annual runoff:
rann = sum(rv,1);
rann = rann(:);

% Fit regression:
X=ones(nyr,1);
X=[X,(1:nyr)'];
B = regress(rann,X);
rfit = X*B;

% extrapolate:
np = 0;
clear rprd
for iyear = 2010:2016
  np = np+1;
  rprd(np,1) = B(1)+B(2)*(iyear-1991);
end;  

YRS=[YR1:YR2];
YPR=[2010:2016];

figure(1); clf;
axes('Position',[0.08 0.4 0.9 0.45]);
plot(YRS,rann);
hold on;
plot(YRS,rfit,'r');
plot(YPR,rprd,'g');

tt{1}='R = a + b*(year-1991)';
tt{2}=sprintf('a=%4.1f, b=%4.1f',B(1),B(2));
text(1992,1050,tt,'Fontsize',11);
title('Annual Total Greenland Runoff, 1992-2010, km3/mo');
set(gca,'ylim',[750 1280],...
	'tickdir','out',...
	'xlim',[1991.5 2016.5],...
	'xtick',[1992:2:2016],...
	'xgrid','on',...
	'ygrid','on');

txtb = 'extrapolate_Greenland_rivers.m';
bottom_text(txtb,'pwd',1);

%keyboard

% Extrapolate from 2010
% and save data
r0_ann = rann(end);

for iyear = 2011:YE
  Rest_ann = B(1)+B(2)*(iyear-1991);
  dR = Rest_ann/r0_ann;
  
  for im = 1:12
    r0 = GR0(im).runoff;
    rEst = r0*dR;

    fprintf('Year=%i, Mo=%i\n',iyear,im);

    if im==1,
      fprintf('Clear GR \n');
      clear GR
      GR = struct;
    end

% Check;
    atot(im,1) = nansum(nansum(rEst));
    
    GR(im).Note='Extrapolated Greenland Runoff, Used Trend';
    GR(im).year=iyear;
    GR(im).units='km3/mo';
    GR(im).runoff=rEst;
    if im==12 & f_save==1;
      rannE = sum(atot);
      dRe = rannE/r0_ann;
      fprintf(' Annual Greenland runoff: %6.2f km3, Estimated: %6.2f\n',...
	      rannE,Rest_ann);
      fout=sprintf('%sGreenland_runoff_monthly-%4.4i.mat',PTH.river,iyear);
      fprintf('Saving runoff: %s\n',fout);
      save(fout,'GR');
    end

    fplt=logical(0);
% Convert km3/mo to m3/sec:
    km3mo=1; % 1km3/mo
    ndays=30.4167;
    mSv=km3mo*1e9/(3600*24*ndays)*1e-3;  % mSv

    if fplt
      A = rEst;
      A(A==0)=nan;
      clf
      pcolor(A); shading flat;
      axis('equal');
      hold on;
      contour(LN,[-100:5:100],'Color',[0.4 0.4 0.4]);
      contour(LT,[40:5:90],'Color',[0.4 0.4 0.4]);
      caxis([0 0.1]);
      set(gca,'color',[0 0 0]);
      set(gcf,'Color',[1 1 1]);
      stt=sprintf('Runoff, km3/mo, %4.4i %2.2i',YC,nm);
      title(stt,'Fontsize',14);
      colorbar;
    end

  end
end;







