% Mean T/S profiles for Eurasian and Candian basins
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

s_fig = 0;
%pfld  = 'temp';
pfld  = 'salin';
fld0  = pfld;
f_save = 2;  % 1 - start from day 1, =2 - load saved

YR1 = 2018;
YR2 = 2020;
regn = 'ARCc0.04';
ires = 0.04;
expt = 23;
expt_nm = 'dltEddTmltEAPJRA';

hg    = 2^100;
rg    = 9806;

pthout = '/nexsan/people/ddmitry/hycom/data_extracted/';


if f_save == 0
  fprintf(' !!!!!!!!!!!!!!!!! \n');
  fprintf(' Output not saved f_save=%i \n\n',f_save);
end

pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';
ftopo = sprintf('%sdepth_ARCc0.04_17DD.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(HH);
% 
% 
YRPLT = [];
cc=0;
for iyr=YR1:YR2
  for idd=1:7:365
    if idd==1, idd=2; end;
    dJ1  = datenum(iyr,1,1);
    cc=cc+1;
    dnmb = dJ1+idd-1;
    YRPLT(cc,1)=iyr;
    YRPLT(cc,2)=idd;
  end
end

fprintf('Extracting %s %i-%i 0.04fHYCOM-CICE expt%3.3i\n',pfld,YR1,YR2,expt);


[mm,nn]=size(HH);
[II,JJ] = meshgrid([1:nn],[1:mm]);
IJEU = [778        1069
        1169        1812
        1582        1200
        1484         625
        1269         709
        1021         989
         833         938];
IJEU = IJEU*2;

IJCA = [778        1069
        1169        1812
         701        1936
         527        1892
         320        1557
         754        1084];
IJCA = IJCA*2;

dsub = 20;
IN = inpolygon(II,JJ,IJEU(:,1),IJEU(:,2));
IN0 = IN*0;
IN0(1:dsub:end,1:dsub:end) = IN(1:dsub:end,1:dsub:end);
INEU = find(IN0==1 & HH<-10.);
IN = inpolygon(II,JJ,IJCA(:,1),IJCA(:,2));
IN0 = IN*0;
IN0(1:dsub:end,1:dsub:end) = IN(1:dsub:end,1:dsub:end);
INCA = find(IN0==1 & HH<-10.);

fmat_out = sprintf('%shycom%3.3i_%3.3i_meanprof_%s.mat',...
                    pthout,ires*100,expt,pfld);

np=size(YRPLT,1);
cnc=0;
ip1=1;
AEUsum = [];
ACAsum = [];
nrc = 0;  
dlast = 0;

if f_save == 2
  fprintf('Start from saved %s\n',fmat_out);
  load(fmat_out);

  TM = PROF.TM;
  dlast = PROF.TM(end);
  nrc = length(TM);
  
  AEUsum = PROF.prof_euras*nrc;
  ACAsum = PROF.prof_canad*nrc;
  fprintf('Last saved record #%i %s\n',nrc, datestr(dlast));

end

for ip=ip1:np
  tic;
  yr=YRPLT(ip,1);
  iday=YRPLT(ip,2);

  if expt==22
    pthbin = sprintf('/nexsan/archive/ARCc0.04_022/data/%4.4i/',yr);
    fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
    finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);
  else
    pthbin = sprintf('/nexsan/people/ddmitry/hycom/ARCc0.04_%3.3i/data/%4.4i_%s/',...
                    expt,yr,expt_nm);
    fina = sprintf('%s%3.3i_archv.%4.4i_%3.3i_00.a',pthbin,expt,yr,iday);
    finb = sprintf('%s%3.3i_archv.%4.4i_%3.3i_00.b',pthbin,expt,yr,iday);
  end

  if ~exist(fina,'file');
    fprintf('Not found: %s\n\n',fina);
    continue;
  end

  cnc=cnc+1;
  dnmb=datenum(yr,1,1)+iday-1;
  DV=datevec(dnmb);

  if dnmb <= dlast
    fprintf('%i %s <= last saved %s, skipping\n',cnc,datestr(dnmb),datestr(dlast));
% Check load data counter:
    if dnmb == dlast & cnc ~= nrc
      error('Check day count in saved mat file')
    end
    continue;
  end

  fprintf(':: %4.4i_%2.2i_%2.2i: %s\n',DV(1:3),fina);

  [F,n,nlev] = read_hycom(fina,finb,pfld);
  F(F>1e6)=nan;

  AEU=squeeze(F(:,INEU));
  [ae1,ae2]=size(AEU);

  ACA=squeeze(F(:,INCA));
  [ac1,ac2]=size(ACA);

%
% Prepare vertical thicknesses 
  fld='thknss';
  [F,n,m,l] = read_hycom(fina,finb,fld);
  F(F>1e10)=nan;
  F(F<0.1)=nan;
  F=F/rg;
  DEU=squeeze(F(:,INEU));
  DEU(DEU==0)=1.0e-16;  % for interpolation, avoid nans
  DCA=squeeze(F(:,INCA));
  DCA(DCA==0)=1.0e-16;  % for interpolation, avoid nans


% Interpolate into reg. vert. grid:
  fprintf('Interpolating into Z-levels ...\n');
  ZZi=[0:-1:-100,-105:-5:-500, -510:-10:-2000, -2020:-20:-6000]';

  for iprf = 1:2
    if iprf == 1
      AA = AEU;
      Dsec = DEU;
    else
      AA = ACA;
      Dsec = DCA;
    end

    clear ZZb
    Dsec(isnan(Dsec))=0;
    ZZb(1,:)=-Dsec(1,:);
    for kk=2:l
      ZZb(kk,:)=ZZb(kk-1,:)-Dsec(kk,:);
    end

    [nl,npb]=size(ZZb);
    ZZ=zeros(nl+1,npb);
    ZZ(2:nl+1,:)=ZZb;

% Depths of middle of the cells:
    ZM = [];
    ZM(1,:)=0.5*ZZ(1,:);
    for kk=1:l
      ZM(kk,:)=0.5*(ZZ(kk+1,:)+ZZ(kk,:));
    end

    npb = size(AA,2);
    nli = length(ZZi);
    AAi = zeros(nli,npb)*nan;

    for ipp=1:npb
      Ib=min(find(Dsec(:,ipp)<=1.e-9));
      if Ib==1, continue; end;
      tt=AA(:,ipp);
      zz=ZM(:,ipp);
      hb=-sum(Dsec(:,ipp));
      ibz = max(find(ZZi>=hb));
      for kl=Ib:nl
        zz(kl)=zz(kl-1)-0.1;
      end;
      zz=[0;zz];
      tt=[tt(1);tt];
      if zz(nl)>ZZi(end)
        zz(nl)=ZZi(end);
      end
      aai = interp1(zz,tt,ZZi,'pchip');
      aai(ibz+1:end)=nan;
      AAi(:,ipp)=aai;
    end

    if iprf == 1
      AEUi = AAi;
    else
      ACAi = AAi;
    end

  end

  if cnc == 1
    AEUsum = AEUi;
    ACAsum = ACAi;
  else
    AEUsum = AEUsum + AEUi;
    ACAsum = ACAsum + ACAi;
  end

  TM(cnc,1) = dnmb;
  PrfEU = nanmean(AEUsum/cnc,2);
  PrfCA = nanmean(ACAsum/cnc,2);

  PROF.TM = TM;
  PROF.prof_euras = PrfEU;
  PROF.prof_canad = PrfCA;
  PROF.Depths = ZZi;

  if mod(cnc,10) == 0 & f_save > 0
    fprintf('Saving output ---> %s\n',fmat_out);
    save(fmat_out,'PROF');
  end
   
  f_plt = 0;
  if f_plt == 1
    figure('Position',[1639         491         864         850]);
    axes('Position',[0.1, 0.1, 0.35, 0.8]);
    hold on
    plot(PrfEU,ZZi)
    set(gca,'xgrid','on',...
            'ygrid','on');
    title(sprintf('Mean %s Euras',pfld));

    axes('Position',[0.55, 0.1, 0.35, 0.8]);
    hold on
    plot(PrfCA,ZZi)
    set(gca,'xgrid','on',...
            'ygrid','on');
    title(sprintf('Mean %s Canad',pfld));

    keyboard
  end 

  fprintf('Min/ max %s: %5.2f/%5.2f, %5.2f/%5.2f\n', ...
          pfld, nanmin(PrfEU), nanmax(PrfEU), nanmin(PrfCA), nanmax(PrfCA));
  fprintf('Processed 1 record %6.2f min\n',toc/60);

end

if f_save > 0
  fprintf('Saving output ---> %s\n',fmat_out);
  save(fmat_out,'PROF');
end

fprintf('  ALL DONE  \n\n');



