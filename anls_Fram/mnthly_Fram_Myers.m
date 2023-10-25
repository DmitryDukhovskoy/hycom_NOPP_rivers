% Calculate monthly means
% form 5-day means
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

YR1=2002;
YR2=2016;

pthdat='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/data_Myers/FRAM/';

btx = 'calc_FWFlux_Fram_Myres.m';

iSp=199; % index of Spitsbergen coast
mold=13;
for iyr=YR1:YR2
  dJ1=datenum(iyr,1,1);

  vmm=[];
  smm=[];
  tmm=[];
  TMv=[];
  TMt=[];
  ccv=0;
  cct=0;

  for iday=5:1:365
    dnmb=dJ1+iday-1;
    dv=datevec(dnmb);
    mo=dv(2);
    mday=dv(3);
    fnm=sprintf('%sANHA12-EXH006_y%4.4im%2.2id%2.2i_FRAMclip_gridV.nc',...
	      pthdat,iyr,mo,mday);
    fprintf('Processing: %s\n',datestr(dnmb));
    
    if ~exist('ZZ','var');
      ZZ=-nc_varget(fnm,'depthv');
      xlon=nc_varget(fnm,'nav_lon');
      xlat=nc_varget(fnm,'nav_lat');
      xlon=xlon(1:iSp);
      xlat=xlat(1:iSp);
      nz=length(ZZ);
      nx=length(xlon);
    end

    if exist(fnm,'file')
      ccv=ccv+1;
      V=squeeze(nc_varget(fnm,'vomecrty'));
      V=V(:,1:iSp);
      vmm(ccv,:,:)=V;
      TMv(ccv)=dnmb;
    end

% T &S    
    fnm=sprintf('%sANHA12-EXH006_y%4.4im%2.2id%2.2i_FRAMclip_gridT.nc',...
	      pthdat,iyr,mo,mday);
    
    if exist(fnm,'file')
      S=squeeze(nc_varget(fnm,'vosaline'));
      S=S(:,1:iSp);
    
      T=squeeze(nc_varget(fnm,'votemper'));
      T=T(:,1:iSp);

      cct=cct+1;
      smm(cct,:,:)=S;
      tmm(cct,:,:)=T;
      TMt(cct)=dnmb;
    end
    
%    pcolor(xlon,ZZ,V); shading flat; colorbar
%    pcolor(xlon,ZZ,T); shading flat; colorbar
%    pcolor(xlon,ZZ,S); shading flat; colorbar
%keyboard    
  end
%  keyboard
% Monthly averages:
  DVt=datevec(TMt);
  DVv=datevec(TMv);
  iic=iyr-YR1+1;
  
  for im=1:12
    I=find(DVv(:,2)==im);
    A=vmm(I,:,:);
    VTS(iic).V(im,:,:)=squeeze(nanmean(A,1));

    I=find(DVt(:,2)==im);
    A=tmm(I,:,:);
    VTS(iic).T(im,:,:)=squeeze(nanmean(A,1));
    A=smm(I,:,:);
    VTS(iic).S(im,:,:)=squeeze(nanmean(A,1));
    if im==1
      VTS(iic).Z=ZZ;
      VTS(iic).X=xlon;
      VTS(iic).Year=iyr;
    end
  end
    
end

fout=sprintf('%smonthly_VTS.mat',pthdat);
fprintf('Saving %s\n',fout);
save(fout,'VTS');
