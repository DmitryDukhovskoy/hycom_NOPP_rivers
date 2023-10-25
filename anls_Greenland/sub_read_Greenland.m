function [dGR,yr] = sub_read_Greenland;
% Read Greenland runoff - Bamber data
% calculate surplus flux
% wrt pre-1990 mean
fprintf('Better to use sub_read_Greenland_v3.m !!!! \n');
fprintf('Reading Greenland runoff fields ...\n');
PTH.river='/Net/ocean/ddmitry/arctic_AOregimes/data/Greenland_rivers/';

YearGr=2004;
mday=[31;28;31;30;31;30;31;31;30;31;30;31];
if mod(YearGr,4)==0,
  mday(2)=29;
end
fgrgr=sprintf('%sGreenland_grid.mat',PTH.river);
load(fgrgr);
LN=GRgrd.LN;
LT=GRgrd.LT;

% Greenland:
rv61_90=876; % mean total runoff over 1961-1990
             % estimate as mean(rv), for YR1=1961, YR2=1990;

YR1=1958;
YR2=2010;
nyr=0;
clear rv Rr
for iyear=YR1:YR2
  YearGr=iyear;
  fgrgr=sprintf('%sGreenland_grid.mat',PTH.river);
  fgrrv=sprintf('%sGreenland_runoff_monthly-%4.4i.mat',PTH.river,YearGr);
  fprintf('===>  Reading %s\n',fgrrv);
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
  end
end;

atot=sum(rv,1); % annual runoff km3/yr
Rmin=min(rv,[],2); % min runoff over N years
Rmax=max(rv,[],2);

%keyboard

YR = [1958:2010];

% Analysis of cumulative river runoff
% using all data: need to change YR1 to 1958 above
yr=[1993:YR2];
a1=atot(1:33); % 1958-1990
a2=atot(36:end); % 1993-2010
cs=sum(a2-mean(a1));
% For 2010-2016 assume 2010
ae=a2(end);
a2=[a2,ae,ae,ae,ae,ae,ae];
yr=[yr,2011,2012,2013,2014,2015,2016];
dGR=a2-mean(a1); % Greenland flux anomaly, km3/yr
  	     

return