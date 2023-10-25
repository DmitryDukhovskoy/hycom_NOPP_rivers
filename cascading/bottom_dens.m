% Calculate monthly mean near-bottom densities
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

yr1=2007;
yr2=2009;

regn = 'ARCc0.08';
%expt = 110; % no Greenland runoff  
expt = 112;  % Greenland runoff

rg = 9806;
% Thershold dRho
dRho0=1e-4*1023;


pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat2/',expt);

s_mat = 2; % =2 - load and start from last saved


ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
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


fprintf('Bottom Density %i-%i\n',yr1,yr2);

Lmsk       = HH*0;
Lmsk(HH<0) = 1;

Iocn=find(Lmsk==1);
nIocn=length(Iocn);
Ilnd=find(Lmsk==0); 

YPLT=[];
cc=0;
dday=7;
for yr=yr1:yr2
  d1=datenum(yr,1,1);
  d2=datenum(yr,12,31);
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
RBTM=struct;
mnxt=-1;
for icc=1:nplt
  yr=YPLT(icc,1);
  mo=YPLT(icc,2);
  mday=YPLT(icc,3);
  dnmb=YPLT(icc,4);
  iday=dnmb-datenum(yr,1,1)+1;
%  if dend==0
%    dend=min([datenum(yr+1,4,1)-dday,YPLT(end,4)]);
%  end
  fmat = sprintf('%sarc08_%3.3i_btmRho_%i.mat',pthmat,expt,yr);

  if mnxt<=0, 
    for im=1:12
      RBTM(im).Rho=HH*0;
      RBTM(im).Cntr=HH*0;
    end
  end   

  if icc<nplt
    mnxt=YPLT(icc+1,2);
  else
    mnxt=mo+100;
  end

  if s_mat==2
    if icc==1 
      load(fmat);
    else
      if YPLT(icc-1,1)<YPLT(icc,1)
        load(fmat);
      end
    end

    rr=RBTM(mo).Rho;
    if max(max(rr))>=999.9
      fprintf('%i/%i processed, skipping ...\n',yr,mo);
      continue
    end
  end
  s_mat=1;

%  pthbin = '/Net/mars/ddmitry/hycom/ARCc0.04/output/';
%  pthbin = sprintf('/nexsan/hycom/ARCc0.04_%3.3i/data/%i/',expt,yr);  
  pthbin = sprintf('/nexsan/archive/%s_%3.3i/data/%i/',regn,expt,yr);
  if expt~=110,
    pthbin = sprintf('/nexsan/hycom/%s_%3.3i/data/%i/',regn,expt,yr);
  end


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
  [F,n,m,vlv] = read_hycom(fina,finb,fld);
  F(F>1e20)=nan;                        
  S=F;  

  fld='temp';
  [F,n,m,l] = read_hycom(fina,finb,fld); 
  F(F>1e20)=nan;                         
  Temp=F; 

%  keyboard

  fprintf('Starting bottom RHO calculation ...\n');
  
  Rho=HH*0;
  cnc=HH*0;
  for k=1:vlv
%    fprintf('level=%i\n',k);
    z1=squeeze(ZZ(k+1,:,:));
    z0=squeeze(ZZ(k,:,:));
    dbt=abs(z1-HH);
    dz=abs(z0-z1);
    Ib=find(dbt<0.1 & dz>0.1);
    if isempty(Ib), continue; end;
    ss=squeeze(S(k,:,:));
    tt=squeeze(Temp(k,:,:));

    rho=sw_dens0(ss,tt);
    Rho(Ib)=rho(Ib);
    cnc(Ib)=cnc(Ib)+1;
  end
  RBTM(mo).Rho=RBTM(mo).Rho+Rho;
  RBTM(mo).Cntr=RBTM(mo).Cntr+cnc;  

  fprintf(' Processing time: %6.2 min\n',toc/60);

% End of month
  if mo~=mnxt
    A=RBTM(mo).Rho;
    cnc=RBTM(mo).Cntr;
    I=find(cnc==0);
    A(I)=nan;
    A=A./cnc;
    RBTM(mo).Rho=A;
    fprintf('mo=%i, Mean btm dns=%6.2f\n',mo,nanmean(nanmean(A)));
    fprintf('    Min btm dns=%6.2f, max =%6.2f\n',...
            min(min(A)),max(max(A)));
    fprintf('Min cntr=%i, max cntr=%i\n',min(min(cnc)), max(max(cnc)));
 
% Save end of year/period
%    fmat = sprintf('%sarc08_%3.3i_btmRho_%i.mat',pthmat,expt,yr);
    fprintf('Saving end-of-month %s\n\n',fmat);
    save(fmat,'RBTM');
    if mo==12 | mnxt>20 
      fprintf('Clearing arrays for next year\n');
      for im=1:12
        RBTM(im).Rho=HH*0;
        RBTM(im).Cntr=HH*0;
      end
    end
  end  % if mo

end



  
  
