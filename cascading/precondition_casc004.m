% Extract T/S at some locations
% and calculate dRho = surf-bottom
% to identify locations of surf-to=bottom mixing
% 
% To be able to be cascaded downslope in the benthic layer, 
% the mixed layer depth (MLD) on the shelf-slope must reach
% the bottom
%
% Monthly mean T/S averaged within the layers
% for 2 experiments (with & without Greenland runoff)
% expt_110 - no Greenland runoff
% expt 112 - with Greenland runoff
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

yr1=2006;
yr2=2007;

regn = 'ARCc0.04';
expt = 011; % no Greenland runoff  
%expt = 012;  % Greenland runoff

rg = 9806;
% Thershold dRho
dRho0=1e-4*1023;


pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
%pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/topo_grid/';
pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/%s/%3.3i/data_mat/',regn,expt);

s_mat = 1; % =2 - load and start from last saved


ftopo = sprintf('%sdepth_ARCc0.04_17DD.nc',pthtopo); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[XX,YY]=meshgrid((1:nn),(1:mm));


fprintf('expt %3.3i\n\n',expt);
if s_mat==0
  fprintf('!!!! Mat file is not saved !!!!! \n');
  fprintf('!!!! Mat file is not saved !!!!! \n');
else
  fprintf('Mat file saved is ON\n');
end


fprintf('Convection Cascading %i-%i\n',yr1,yr2);

% Check shelf regions only, shallower than 500 m
% Mask of the region of interest
% exclude Pacific Ocean and North.Atl.
IJ=[632        1927
         546        1890
         586        1773
         391        1664
         292        1475
         231        1189
         261         984
         280         563
         372         298
         442         148
         536         153
         474         394
         472         558
         589         366
         858         361
        1005         348
        1210         358
        1265         527
        1350         683
        1457         610
        1561         703
        1583        1161
        1567        1303
        1375        1573
        1226        1769
        1081        1856
         836        1940
         663        1921];

IJ=IJ*2;

INP=inpolygon(XX,YY,IJ(:,1),IJ(:,2));
Lmsk     = HH*0;
Lmsk(HH>-800 & HH<0 & INP==1) = 1;

Iocn=find(Lmsk==1);
nIocn=length(Iocn);
Ilnd=find(Lmsk==0); 

% Winter only Oct Year1 - March 31 Year1+1
YPLT=[];
cc=0;
dday=5;
for yr=yr1:yr2
  d1=datenum(yr,10,1);
  d2=datenum(yr+1,3,31);
  if yr+1>2007   % last year of the experiment
    d2=datenum(2007,12,31);
  end
  
  for dd=d1:dday:d2
    cc=cc+1;
    DV=datevec(dd);
    YPLT(cc,1)=DV(1);
    YPLT(cc,2)=DV(2);
    YPLT(cc,3)=DV(3);
    YPLT(cc,4)=dd;
  end
end
nplt=cc;


cnc=0;
CONV=HH*0;
Tconv=HH*0;
Sconv=HH*0;
MXD=struct;
MXD.Conv_stat=CONV;
MXD.Tconv=Tconv;
MXD.Sconv=Sconv;
dend=0;
for icc=1:nplt
  yr=YPLT(icc,1);
  mo=YPLT(icc,2);
  mday=YPLT(icc,3);
  dnmb=YPLT(icc,4);
  iday=dnmb-datenum(yr,1,1)+1;
  if dend==0
    dend=min([datenum(yr+1,4,1)-dday,YPLT(end,4)]);
  end
  


%  pthbin = '/Net/mars/ddmitry/hycom/ARCc0.04/output/';
%  pthbin = sprintf('/nexsan/hycom/ARCc0.04_%3.3i/data/%i/',expt,yr);  
%  pthbin = sprintf('/nexsan/archive/%s_%3.3i/data/%i/',regn,expt,yr);
  pthbin = sprintf('/nexsan/hycom/ARCc0.04_%3.3i/data/%4.4i/',expt,yr);

  fprintf('\n %i/%i/%i\n',yr,mo,mday);
  fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
  finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);

  if ~exist(fina,'file');
    fprintf('Not found: %s\n\n',fina);
    continue;
  end

  tic; 

  [ZM,ZZ] = sub_zz_zm(fina,finb,HH);

  fld='salin';
  [F,n,m,l] = read_hycom(fina,finb,fld);
  F(F>1e20)=nan;                        
  S=F;  

  fld='temp';
  [F,n,m,l] = read_hycom(fina,finb,fld); 
  F(F>1e20)=nan;                         
  Temp=F; 

%  keyboard

  fprintf('Starting MLD calculation ...\n');

  ibb=0;
  mbb=0;
  for ii=1:nIocn   
    if mod(ii,50000)==0
      rr=ii/nIocn*100;
      fprintf('Done %3.1f%% ...\n',rr);
    end;
    id=Iocn(ii);
    [j,i]=ind2sub(size(HH),id);

    Tz=squeeze(Temp(:,id));
    Sz=squeeze(S(:,id));
    ZMz=squeeze(ZM(:,id));
    Zz=squeeze(ZZ(:,id));
    dZ=abs(diff(Zz));
    Ibb=find(dZ<1e-12); % bottom - can creates NaN in the middle
    di=max(diff(Ibb));
    if di>1
      for iB=1:length(Ibb)-1
	di=Ibb(iB+1)-Ibb(iB);
	if di>1
	  dZ(Ibb(iB))=1e-6;
	end
      end
      Ibb=find(dZ<1e-12); % bottom - can creates NaN in the middle
    end
    Tz(Ibb)=nan;
    Sz(Ibb)=nan;
% Check bottom
    cdz=round(cumsum(dZ)*1e4)/1e4;
    Hb=round(abs(HH(id))*1e4)/1e4;
    dhb=abs(Hb-cdz);
    ibb=min(find(dhb==min(dhb)));
    if abs(abs(HH(id)/cdz(ibb))-1)>1e-2
      fprintf('Bottom mismatch: %8.4f vs %8.4f ...\n',abs(HH(id)),cdz(ibb));
      keyboard
    end

    Tz(ibb:end)=nan;
    Sz(ibb:end)=nan;
    if isnan(Tz(1)) % land, possible mismatch btw my topo and Global Reanalysis
      continue;
    end

    [Zmld,dR] = sub_Dlt_Rho(Tz,Sz,ZMz,Zz,ibb,Hb,dRho0);

    if abs(abs(Zmld)-abs(Hb))/abs(Hb)<0.05, 
%      fprintf('Convection point: %6.1f i=%i, j=%i\n',Hb,i,j);
      CONV(id)=CONV(id)+1;
      Tconv(id)=Tconv(id)+nanmean(Tz(1:ibb-1));
      Sconv(id)=Sconv(id)+nanmean(Sz(1:ibb-1));
      mbb=mbb+1;
    end

%keyboard
  end
  cnc=cnc+1;
  fprintf('Recrds %i, Found Convection points: %i\n',cnc,mbb);
  fprintf('====  All ocean pnts: %8.1f min\n',toc/60);
  
  MXD.Info='Surf-to-bottom convection: # per winter, aver T, S';

  MXD.Conv_stat=CONV;
  MXD.Tconv=Tconv;
  MXD.Sconv=Sconv;
  MXD.Nrecords=cnc;
  MXD.TM(cnc)=dnmb;

  f_plt=0;
  if f_plt==1
    figure(10); clf;
    hold on;
    contour(HH,[0 0],'k');
    contour(HH,[-500 -500],'b');
    Ic=find(CONV==1);
    plot(XX(Ic),YY(Ic),'g.');
  end

% Save at the end of winter
  if dnmb>=dend
    fprintf(' End of Winter %s ...\n',datestr(dnmb));
    dv=datevec(dnmb);
    if dv(2)<10,
      yrM=dv(1)-1;
    else
      yrM=dv(1);
    end
    
    fmat = sprintf('%sarc04_%3.3i_precond_wint%i.mat',pthmat,expt,yrM);
    
    TT=MXD.Tconv;
    SS=MXD.Sconv;
    CONV=MXD.Conv_stat;
    II=find(CONV>0);
    TT(II)=TT(II)./CONV(II);
    SS(II)=SS(II)./CONV(II);
    Inan=find(CONV==0);
    TT(Inan)=nan;
    SS(Inan)=nan;
    MXD.Tconv=TT;
    MXD.Sconv=SS;
  
    I0=find(CONV==0);
    TT(I0)=nan;
    SS(I0)=nan;
    tm=nanmean(nanmean(TT));
    sm=nanmean(nanmean(SS));
    tmin=min(TT(II));
    tmx=max(TT(II));
    smin=min(SS(II));
    smx=max(SS(II));
    fprintf('Mean conv T=%6.2f, min/max T=%6.2f / %6.2f\n',tm,tmin,tmx); 
    fprintf('Mean conv S=%6.2f, min/max S=%6.2f / %6.2f\n',sm,smin,smx); 
    fprintf('Mean Prob. convec=%5.4f\n',mean(CONV(II))/cnc);
    
    if s_mat==1
      fprintf(':::: ==> Saving %s\n',fmat);
      save(fmat,'MXD');
    end
  
    fprintf(' Zeroing fields ...\n');
    CONV=HH*0;
    Tconv=HH*0;
    Sconv=HH*0;
    cnc=0;
    MXD.Conv_stat=CONV;
    MXD.Tconv=Tconv;
    MXD.Sconv=Sconv;
    dend=min([datenum(yr+1,4,1)-dday,YPLT(end,4)]);
    
    fprintf(' ==========  New cycle Ends: %s\n',datestr(dend));
    
  end
  
end

%if s_mat==1
%  fprintf('End of run: Saving %s\n',fmat);
%  save(fmat,'CONV');
%end


  
  