addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

Sref=34.9;
YR1=2002;
YR2=2016;


pthdat='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/data_Myers/FRAM/';

btx = 'FW_vol_fluxes_myers.m';

fprintf('NEMO FW and Vol fluxes, Sref=%4.2f, %i-%i\n',Sref,YR1,YR2);


% Fraction of transport not measured:
SBE = Fram_moorAWI;
nf=length(SBE);

for ik=1:nf
  Ymr(ik)=SBE(ik).lat_lon(1);
  Xmr(ik)=SBE(ik).lat_lon(2);
end
xSBE=min(Xmr);



iSp=199; % index of Spitsbergen coast
mold=13;
for iyr=YR1:YR2
  dJ1=datenum(iyr,1,1);

  vfm=[];
  vmm=[];
  TMv=[];
  volD=[];
  fwfD=[];
  ccv=0;

  for iday=5:1:365
    dnmb=dJ1+iday-1;
    dv=datevec(dnmb);
    mo=dv(2);
    mday=dv(3);
    fnm=sprintf('%sANHA12-EXH006_y%4.4im%2.2id%2.2i_FRAMclip_gridV.nc',...
	      pthdat,iyr,mo,mday);
    fsnm=sprintf('%sANHA12-EXH006_y%4.4im%2.2id%2.2i_FRAMclip_gridT.nc',...
	      pthdat,iyr,mo,mday);

    if ~exist(fnm,'file') | ~exist(fsnm,'file'); continue; end;
    
    fprintf('Processing: %s\n',datestr(dnmb));
    
    if ~exist('ZZ','var');
      Z=-nc_varget(fnm,'depthv');
      nzz=length(Z);
      ZM=Z;
      dZ=abs(diff(Z));
      dZ(nzz)=dZ(nzz-1);
      ZZ(1)=0;
      for iz=1:nzz
	ZZ(iz+1)=ZM(iz)-0.5*dZ(iz);
      end

      xlon=nc_varget(fnm,'nav_lon');
      xlat=nc_varget(fnm,'nav_lat');
      xlon=xlon(1:iSp);
      xlat=xlat(1:iSp);
      nz=length(ZZ);
      nx=length(xlon);
      
% Grid sizes:
      for ik=1:nx-1
	x0=xlon(ik);
	y0=xlat(ik);
	x1=xlon(ik+1);
	y1=xlat(ik+1);
	dx=distance_spheric_coord(y0,x0,y1,x1);
	DX(ik)=dx;
      end
      DX(nx)=DX(nx-1);
      
      [DXG,DZG]=meshgrid(DX,dZ);
      Acell=DXG.*DZG;
      
% Segments of mooring obs. locations:
      ixW=max(find(xlon<=xSBE)); % Gr Shelf west of F17 mooring
      x=max(Xmr);
      ixE=min(find(xlon>=x)); % easternmost location in x-sgements arrays only
      
    end

    ccv=ccv+1;
%      VF=squeeze(nc_varget(fnm,'vs'));
%      VF=VF(:,1:iSp);

    V=squeeze(nc_varget(fnm,'vomecrty'));
    V=V(:,1:iSp);
    vmm(ccv,1)=nansum(nansum(V.*Acell)); % m3/sec - vol flux
    volD(ccv,:)=nansum(V.*Acell); % transport 1D along section

% Saved is salt mass flux: S*V      
% To recover FW flux need:
% Vfw=U*(S0-S)/S0=U-S/S0*U
%
%      vfw=V-VF/Sref;
%      vfm(ccv,1)=nansum(nansum(vfw.*Acell));  % m3/s of freshwater 
    S=squeeze(nc_varget(fsnm,'vosaline'));
    S=S(:,1:iSp);
    S(S>Sref)=nan;
    dmm=V.*(Sref-S)./Sref;
    vfm(ccv,1)=nansum(nansum(dmm.*Acell));   % m3/s of freshwater 
    fwfD(ccv,:)=nansum(dmm.*Acell);   % FW flux 1D along section
    fprintf('FWFlux=%8.2f mSv\n',vfm(ccv)*1e-3);

    TMv(ccv)=dnmb;

%    pcolor(xlon,Z,VF); shading flat; colorbar

%keyboard
  end
  
  DVv=datevec(TMv);
  iic=iyr-YR1+1;
  
  FWV.Time(iic)=iyr;
  VF1D(iic).Year=iyr;
  VF1D(iic).ixW=ixW;
  VF1D(iic).ixE=ixE;
  for im=1:12
    I=find(DVv(:,2)==im);
    A=vfm(I);
    FWV.FW_flx(iic,im)=nanmean(A,1);
    
    A=vmm(I);
    FWV.Vol_flx(iic,im)=nanmean(A,1);

    VF1D(iic).FWF_1D(im,:)=nanmean(fwfD(I,:));
    VF1D(iic).VolF_1D(im,:)=nanmean(volD(I,:));
  end
  
end


fout=sprintf('%smonthly_FW_vol_fluxes.mat',pthdat);
fprintf('Saving %s\n',fout);
save(fout,'FWV','VF1D');

