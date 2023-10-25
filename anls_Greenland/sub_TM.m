function TD = sub_TM(TM,CYC);
% function TD = sub_TM(TM,CYC);
% Create a continuous time series 
% from 2-cycles time series
%  Returns # of days from day 1 of sim=dj1
% 
dj1 = datenum(2004,1,1); % start of the simulation
iy1 = max(find(CYC<2));
if isempty(iy1)
  fprintf('Only expt 061 is saved\n');
  iy1=0;
  d2=datenum(2011,1,1)-datenum(2004,1,1);
  TD=TM-datenum(2004,1,1)+d2;
  return;
end

dv2 = datevec(TM(iy1));
dv1 = datevec(TM(1));
nyrs= dv2(1)-dv1(1)+1;
[a1,a2]=size(TM);
if a1<a2
  TM=TM';
end;
tm1 = TM(1:iy1);
ic=find(CYC>1);
if isempty(ic),
  fprintf('Only 1 experiment is saved: 060 \n');
  TD=tm1;
else
  tm2 = TM(iy1+1:end);
  d2  = TM(iy1)-TM(1); % day offset for cyc 2
  tm1 = tm1-tm1(1); % days since day1 of simulation, cycle 1
  tm2 = tm2-tm2(1); % days since day1 of sim, cyc=2
  TD  = [tm1;tm2+d2+1];
end


return