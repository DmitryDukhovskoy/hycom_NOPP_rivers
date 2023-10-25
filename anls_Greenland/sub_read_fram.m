function FWT = sub_read_fram(shelfM);
% Data from Norwegian Polar Inst.
% de Steur et al., 2018
% Note different moorings for 1997-2002
% and 2003-2016
% different estimates of FW transports
% Also there are 2 time series for 2003-2016
% with and without shelf mooring
% with the shelf mooring - higher FW flux
% shelfM=0 - discard shelf mooring (note much lower FW transport)
% shelfM=1 - with shelf mooring but shorter T Ser. 
%
% Freshwater transport (relative to a reference salinity of 34.9) and 
% gridded data fields of southward velocity and salinity from the East Greenland Current in Fram Strait. 
% The time series are collected by moorings F11, F12, F13, F14, F17 
% from the Fram Strait Arctic Outflow Observatory (Norwegian Polar Institute) 
% and by moorings F9 and F10 (Alfred Wegener Institute). 
% The moorings were at 79N from Sept. 1997 to Aug. 2002 
% and at 78°50'N from Sept. 2002 to Aug. 2015. 
% Only in Sept. 2003 the extra mooring (F17) was added on the shelf.
%
fprintf('========  READING NPI Fram Freshwater Trt ===\n');
fprintf(' Ignore warnings - missing data in 1997\n');

if shelfM==1
  fprintf('   Shelf Mooring time series is selected, shorter T. Ser\n');
  fprintf('   =======   \n');
else
  fprintf(' Note shelf mooring is not selected - low bias in FWT\n');
end

pthdat='/Net/data2/ddmitry/Fram_FWT/';
TM1=[];
TM2=[];
if shelfM==1
  fl2='Fram_Strait_freshwater_transport_2003-2015-v1.0.nc'; % with shelf mooring
  fp2=sprintf('%s%s',pthdat,fl2);
  fwt2=nc_varget(fp2,'FWT'); % mSv
  tm=nc_varget(fp2,'TIME');
  TM=datenum(1950,01,01)+tm-1;
  TM=round(TM);
  fwt=fwt2; %mSv 
else
% Obs no shelf
% 1997-2002 at 79N
  fl1='Fram_Strait_freshwater_transport_1997-2002-v1.0.nc';
% no shelf 2002-2015 at 78.5N
  fl2='Fram_Strait_freshwater_transport_2002-2015-v1.0.nc'; % no shelf mooring
  fp1=sprintf('%s%s',pthdat,fl1);
  fp2=sprintf('%s%s',pthdat,fl2);

  fwt1=nc_varget(fp1,'FWT'); % mSv
  tm=nc_varget(fp1,'TIME');
  TM1=datenum(1950,01,01)+tm-1;
  fwt2=nc_varget(fp2,'FWT'); % mSv
  tm=nc_varget(fp2,'TIME');
  TM2=datenum(1950,01,01)+tm-1;

  TM=[TM1;TM2];
  TM=round(TM);
  TM1=round(TM1);
  TM2=round(TM2);
  fwt=[fwt1;fwt2]; %mSv 
end

DV=datevec(TM);
dmm=diff(DV(:,2));
yr1=DV(1,1);
yr2=DV(end,1);
iday=round(mean(DV(:,3))); % assuming same date every month
%keyboard

clear fwty fwtm
cc=0;
ccy=0;
for iyr=yr1:yr2
  msum=0;
  ccy=ccy+1;
  kk=0;
  for im=1:12
    dnmb=datenum(iyr,im,iday);
    dTM=abs(TM-dnmb);
    I=find(dTM<=1);
    dd1=datenum(iyr,im,1);
    dnxt=datevec(dnmb+35);
    mnxt=dnxt(2);
    ynxt=dnxt(1);
    dd2=datenum(ynxt,mnxt,1);
    mdays=dd2-dd1;
    
    cc=cc+1;
    TMd(cc,1)=dnmb;
    
    
    if isempty(I) | isnan(fwt(I));
      fwtm(cc,1)=nan;
    else
      fwtm(cc,1)=fwt(I)*1e3*3600*24*mdays*1e-9; %mSv->km3/mo
      kk=kk+mdays;
      msum=msum+fwtm(cc,1);
    end
  end
  if kk<180
    fwty(ccy,1)=nan;
  else
    fwty(ccy,1)=msum/kk*365; % annual mean FWT km3/day->km3/yr
  end
end

%keyboard

ymo=[yr1:1/12:yr2+0.99];
YR=[yr1:yr2];
f_chck=0;
if f_chck==1
  figure(1); clf
  subplot(2,1,1);
  hold
  plot(ymo,fwtm,'Color',[0.7 0.7 0.7]);
  for nn=1:length(YR);
    iyr=YR(nn);
    plot([iyr iyr+0.99],[fwty(nn)/12 fwty(nn)/12],'r-');
  end
  set(gca,'xtick',[yr1:yr2],...
	  'xgrid','on');
  title('Fram Freshwater Trt, Sref=34.9, km3/mo');
  
  subplot(2,1,2);
  plot(YR,fwty);
  title('Annual mean FWT Fram, Sref=34.9, km3/yr');
  
end

FWT.Monthly_km3_mo=fwtm;
FWT.Annual_km3_yr=fwty;
FWT.TM=TMd;
FWT.TM_1st_obs=TM1;
FWT.TM_2nd_obs=TM2;
FWT.Years=YR;

    
return