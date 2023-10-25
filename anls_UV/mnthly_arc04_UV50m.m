% Extract and calculate Monthly mean U,V fields
% high-res simulation 0.04
% averaged over exact depth levels zz1:zz2
%
% see anls_Greenland/monthly_mean_flds.m
% exact averaging: extr_MassTrcr_month.m
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

YR1 = 2007;
YR2 = 2007;

s_mat = 1; % =0 - do not save mat file
           % =1 - save tracer and overwrite existing mat
	   % =2 - skip months where mat file exist

if s_mat==0,
  fprintf('Mat file is not created\n');
elseif s_mat == 1
  fprintf('Mat file will be saved, old mat file will be overridden\n');
elseif s_mat == 2
  fprintf('Extraction is skipped for months where old mat files exist\n');
end

% Averaged over the depth layers:
zz1=0;
zz2=-50;

regn = 'ARCc0.04';
expt = 011;
rg = 9806;


pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';
pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/%s/%3.3i/data_mat2/',regn,expt);

fprintf('%s Monthly mean UV %i-%i in the layer: %i- %im\n',regn,YR1,YR2,zz1,zz2);

meanUV=struct('U',0,'V',0,'cnn',0,'TM',0);
meanUV(1).Fields='monthly mean U in 1 layer';
meanUV(1).Units='m/s';

for kk=2:12
  meanUV(kk)=meanUV(1);
end;

rg=9806;  % convert pressure to depth, m

%figure(1); clf;
%set(gcf,'Visible','off');

%YPLT=[2004,299;2005,351;2006,295;2007,354;2008,286;2009,236];
% Animation:
YRPLT=[];
cc=0;
dday = 7;
for iyr=YR1:YR2
  for idd=1:dday:365
%    if idd==1, idd=2; end;
    cc=cc+1;
    YRPLT(cc,1)=iyr;
    YRPLT(cc,2)=idd;
  end
end

np=size(YRPLT,1);

%ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
ftopo = sprintf('%s/depth_%s_17DD.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

hmsk=HH;
hmsk(HH<0)=nan;

if ~exist('DX','var')
  [DX,DY]=sub_dx_dy(LON,LAT);
end

% Get HYCOM fields:
ip1=1;
mold = 0;

if s_mat==2
  for yr=YR1:YR2
    fmat = sprintf('%smnthUV_%4.4i-%4.4i_%i.mat',pthmat,abs(zz1),abs(zz2),yr);

    fprintf('Loading saved %s\n',fmat);
    load(fmat);
    TMl=1-dday;
    fall = 1;
    for imL=1:12
      dmm = meanUV(imL).U;
      TM  = meanUV(imL).TM;
      if min(min(dmm))==0 & max(max(dmm))==0, 
	fall = 0;
	break; 
      end;
      if TM==0, TM=datenum(yr,imL,28); end;
      TMl = TM(end);
    end
    if fall==1, 
      fprintf('Year %i is complete, skipping to next year...\n',yr);
      continue; 
    end;
    dv=datevec(TMl);
    TMnxt = TMl+dday;
    iday = TMnxt - datenum(dv(1),1,1)+1;
    ip1 = max(find(YRPLT(:,2)<iday))+1;
    dnmb=datenum(yr,1,1)+iday-1;
    DV=datevec(dnmb);
    imo=DV(2);
    mold = imo;

    fprintf('::::  Last saved record = %i %s\n',imL-1,datestr(TMl));
    fprintf('::::  Will start from iday=%i, %s, mold=%i\n',...
	    iday,datestr(TMnxt),mold);
    fprintf('::::  YRPLT=%i %i\n',YRPLT(ip1,:));
%    keyboard
  end
end


for ip=ip1:np
  yr=YRPLT(ip,1);
  iday=YRPLT(ip,2);
  fmat = sprintf('%smnthUV_%4.4i-%4.4i_%i.mat',pthmat,abs(zz1),abs(zz2),yr);
  
  dnmb=datenum(yr,1,1)+iday-1;
  DV=datevec(dnmb);
  imo=DV(2);
  fprintf('Processing Month %i\n',imo);
  
  pthbin = sprintf('/nexsan/hycom/%s_%3.3i/data/%i/',regn,expt,yr);

% Save finished month and zero in arrays
  if mold~=imo
    if exist('meanUV','var') & s_mat>0 & mold>0
      cnn = meanUV(mold).cnn;
      if cnn == 0,
	error('# of saved records = 0');
      else
	fprintf('Saving Monthly mean, # of av. rcrds=%i\n',cnn);
      end

      fprintf('%i Saving meanUV month %i ...\n',yr,mold);
%keyboard
      fprintf('Saving %s\n',fmat);
      save(fmat,'meanUV');
    end

    if mold ==12
      for im = 1:12
        A=meanUV(mold).U*0;
        meanUV(im).cnn=0;
	meanUV(im).TM=0;
	meanUV(im).Uav=A;
	meanUV(im).Vav=A;
      end
    end
    
    mold = imo;
    nrec = 0; 
  end


  fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
  finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);

  if ~exist(fina,'file');
    fprintf('Not found: %s\n\n',fina);
    continue;
  end
  
  fprintf('%4.4i_%2.2i_%2.2i: %s\n',DV(1:3),fina);

  tic;
  [F,n,m,l] = read_hycom(fina,finb,'u-vel.');
  F(F>1e10)=nan;
  U=squeeze(F);
  
  [F,n,m,l] = read_hycom(fina,finb,'v-vel.');
  F(F>1e10)=nan;
  V=squeeze(F);
  
  fprintf('Getting layer depths ...\n');
  [ZM,ZZ] = sub_zz_zm(fina, finb,HH,'f_btm',1);
  
  dz=abs(zz2-zz1);
  zbtm=-9000;
  Uav = sub_zLayer_average(HH,ZZ,U,zbtm,zz1,zz2); % depth-average field 
  Vav = sub_zLayer_average(HH,ZZ,V,zbtm,zz1,zz2); % depth-average field 
  
  %
% Need to collocate U, V and dH
%  [ll,mm,nn]=size(U);
  Tx=Uav*0;
  Ty=Vav*0;
  for ii=2:nn-1
    u1 = Uav(:,ii);
    u2 = Uav(:,ii+1);
    uu = 0.5*(u1+u2);
%    dh=dH(:,ii);
%    Tx(:,ii)=uu.*dh;
    Tx(:,ii) = uu; % collocated U in the layer
  end

  for jj=2:mm-1
    v1 = Vav(jj,:);
    v2 = Vav(jj+1,:);
    vv = 0.5*(v1+v2);
%    dh = squeeze(dH(jj,:));
%    Ty(jj,:)=vv.*dh;
    Ty(jj,:)=vv;
  end
  
  mo=DV(2);
  A=meanUV(mo).U;
  B=meanUV(mo).V;
  cnn=meanUV(mo).cnn;
  Aav=cnn/(cnn+1)*A+1/(cnn+1)*Tx;
  Bav=cnn/(cnn+1)*B+1/(cnn+1)*Ty;
  cnn=cnn+1;
  meanUV(mo).U=Aav;
  meanUV(mo).V=Bav;
  meanUV(mo).cnn=cnn;
  meanUV(mo).TM(cnn)=dnmb;

  
  s=sqrt(Aav.^2+Bav.^2);
  fprintf(' == Speed: max=%3.2f min=%3.2f\n',max(max(s)),min(s(s>0.01)));
  fprintf('1 record processed in %8.6f min\n\n',toc/60);
  
%  keyboard
  f_pltU = 0;
  if f_pltU>0 & mod(ip,10)==0
    nf=1;
    stl=sprintf('U %4.4i/%2.2i/%2.2i, Layer %i-%i',DV(1:3),abs(zz1),abs(zz2));
    u=Aav(:,:);
    v=Bav(:,:);
    s=sqrt(u.^2+v.^2);
    c1=0;
    c2=0.5;
    xl1=200;
    xl2=1200;
    yl1=100;
    yl2=1300;
    sub_plot_divu(s,LON,LAT,nf,HH,u,v,stl,xl1,xl2,yl1,yl2,c1,c2);
  end

 
end;  % time loop

if s_mat>0
  fprintf('===End:   Saving meanUV month %i ...\n',mold);

%  fmat = sprintf('%smnthUV_lr%2.2i_%i.mat',pthmat,plr,yr);
  fprintf('Saving %s\n',fmat);
  save(fmat,'meanUV','-v7.3');
end




