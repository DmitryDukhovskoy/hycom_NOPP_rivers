% Calculate monthly wind stress curl over the central Arctic
% and wind over the Fram and Nordic Seas
% from NCEP CFSR and CFSv2
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

pthdat1 = '/Net/kronos/ddmitry/ncep_cfsr/';
pthdat2 = '/Net/kronos/ddmitry/ncep_cfsv2/';
pthmat  = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.72/data_mat/';
pth72   = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.72/topo_grid/';

btx = 'calc_Nordic_wndcrl.m';

% Get indices for interpolating onto HYCOM ARCc0.72 grid:
fmat = sprintf('%scfsr_gridindx_arc072_nghb.mat',pthmat);
load(fmat);
mc = INDX.dim_rows;
nc = INDX.dim_colms;

% Topo for ARCc0.72
fsv=[pth72,'new_bath072.mat'];
load(fsv);
LN=elon;
LT=alat;
HH=hnew;
HH(isnan(HH))=10;
clear hnew elon alat;
[mm,nn]=size(HH);
[DX,DY]=sub_dx_dy(LN,LT);

% Get PANG - angles to rotate vectors into HYCOM grid
% b) Rotate the velocities HYCOM ARCc to x-wards y-wards.  
%The array pang in regional.grid can be used to do this.  
%The eastwards,northwards to x-wards,y-wards code is in ALL/force/src/wi.f:
%
%              COSPANG  = COS(PANG(I,J))
%              SINPANG  = SIN(PANG(I,J))
%              TXM(I,J) = COSPANG*TXMIJ + SINPANG*TYMIJ
%              TYM(I,J) = COSPANG*TYMIJ - SINPANG*TXMIJ
fina = sprintf('%sregional.grid.a',pth72);
finb = sprintf('%sregional.grid.b',pth72);
[PANG,nn,mm] = read_pang(fina,finb);
COSPANG = cos(PANG);
SINPANG = sin(PANG);

cc=0;
CRL=[];
for yr=1993:2016
  
  if yr>=2011
    fmat = sprintf('%scfs_gridindx_arc072_nghb.mat',pthmat);
    load(fmat);
    mc = INDX.dim_rows;
    nc = INDX.dim_colms;
  end
  
  if yr<2011
    fnm = sprintf('%scfsr-sea_%i_mon_uv-10m.nc',pthdat1,yr);
  else
    fnm = sprintf('%scfsv2-sec2_%i_mon_uv-10m.nc',pthdat2,yr);
  end

  if ~exist('X','var')
    X  = nc_varget(fnm,'Longitude');
    Y  = nc_varget(fnm,'Latitude');
  end
  dm = nc_varget(fnm,'MT');
  TM=datenum(1900,12,31)+dm;

% Northward and eastward wind components
  U = nc_varget(fnm,'wndewd');
  V = nc_varget(fnm,'wndnwd');
  
  for im=1:12
    dnmb = TM(im);
    DtV = datevec(dnmb);
    fprintf('%i/%2.2i/%2.2i\n',DtV(1:3));
    
    uu=squeeze(U(im,:,:));
    vv=squeeze(V(im,:,:));
    s=sqrt(uu.^2+vv.^2);

    uh = uu(INDX.II);
    vh = vv(INDX.II);
    ur = COSPANG.*uh+SINPANG.*vh;
    vr = COSPANG.*vh-SINPANG.*uh;
    sr = sqrt(ur.^2+vr.^2);
    
% Mean curl over Nordic Seas getpts
    sC=0;
    sA=0;
    ik1=102;
    ik2=122;
    jk1=62;
    jk2=92;
    di=5;
    for ik=ik1:5:ik2
      for jk=jk1:5:jk2
	A=HH(jk-di:jk+di,ik-di:ik+di);
	if max(max(A))>0, continue; end; % only ocean points
	u1=ur(jk-di,ik);
	u2=ur(jk+di,ik);
	v1=vr(jk,ik-di);
	v2=vr(jk,ik+di);
	dx=abs(sum(DX(jk,ik-di:ik+di)));
	dy=abs(sum(DY(jk-di:jk+di,ik)));
	crl = ((v2-v1)/dx-(u2-u1)/dy);
        acell=dx*dy;
	sC=sC+crl*acell;
	sA=sA+acell;
      end
    end
    crl_mean = sC/sA; % s-1
    cc=cc+1;
    CRL(cc,1)=crl_mean;
    DNM(cc,1)=dnmb;
    

% Mean North/South wind through Fram Strait
    jF=104;
    iF1=101;
    iF2=118;
    vmm=vr(jF,iF1:iF2);
    VFr(cc,1) = mean(vmm);
    
    s_chck=0;
    if s_chck>0
      sub_plotWnd_chck(sr,ur,vr,nn,mm,s,uu,vv,...
		       LN,LT,HH,INDX,DtV,mc,nc);      
      keyboard
    end
    
    
  end
  
end

% ==================
% Plot Beaufort Gyre
% ==================
DV=datevec(DNM);
nc=length(CRL);
TT=[0:nc-1]'/12+DV(1);

% Winter vort:
IS=find(DV(:,2)>4 & DV(:,2)<10);
IW=find(DV(:,2)<=4 | DV(:,2)>9);

cp=0;
for iyr=1993:2016
% Winter
%  IYW=find(DV(:,1)==iyr & (DV(:,2)<4 | DV(:,2)>9));
  IYW=find((DV(:,1)==iyr-1 & DV(:,2)>9) |...
	   (DV(:,1)==iyr & DV(:,2)<4));
  dmm=mean(CRL(IYW));
  vmm=mean(VFr(IYW));
  cp=cp+1;
  CW(cp,1)=dmm; % winter curl by years
  TY(cp,1)=iyr;
  VFrW(cp,1)=vmm;
  
% Summer
  IYS=find(DV(:,1)==iyr & (DV(:,2)>4 & DV(:,2)<10));
  smm=mean(VFr(IYS));
  VFrS(cp,1)=smm;
end



irr=find(TT<1997);
ipp=find(TT>=1997);
mn90=mean(CRL(irr)); % annual mean curl prior 1997
mn00=mean(CRL(ipp)); % annueal mean curl after 1997

iw1=find(TY<1997);
iw2=find(TY>=1997);
wmn90=mean(CW(iw1)); % winter mean curl prior 1997
wmn00=mean(CW(iw2)); % winter mean curl >1997

%For it=1:nc
figure(3); clf;
axes('Position',[0.08 0.58 0.85 0.35]);
hold on;
yl1 = -8e-6;
yl2 = 15e-6;
for iyr=1993:2017
  ix1=iyr-0.25;
  ix2=iyr+0.25;
  xx=[ix1 ix2 ix2 ix1];
  yy=[yl1 yl1 yl2 yl2];
  pp = patch(xx,yy,[0.9 0.99 1]);
  set(pp,'Edgecolor','none');
end
plot(TT,CRL,'k');
plot([TT(1) TT(max(irr))],[mn90 mn90],'r--');
plot([TT(min(ipp)) TT(max(ipp))],[mn00 mn00],'g--');

set(gca,'tickdir','out',...
	'xlim',[TT(1) TT(end)],...
	'xtick',[floor(TT(1)):1:ceil(TT(end))],...
	'ylim',[yl1 yl2],...
	'xgrid','on');
title('Wind curl, Nordic Seas');



axes('Position',[0.08 0.08 0.85 0.35]);
bb=bar(TY,CW);
hold on;
plot([1992.5 1996.25],[mn90 mn90],'r--');
plot([1996.7 2016.5],[mn00 mn00],'g--');

set(bb,'FaceColor',[0. 0.6 0.8],'EdgeColor','none');
set(gca,'tickdir','out',...
	'xlim',[1992.5 2016.5],...
	'xtick',[1993:2016]);
title('Winter Wind curl, Nordic Seas');

bottom_text(btx,'pwd',1);


XV = [-10:2:16]*1e-6;
Hw = hist(CRL(IW),XV);
Hs = hist(CRL(IS),XV);

figure(4); clf;
subplot(2,2,1);
bw=bar(XV,Hw,1);
set(bw,'FaceColor',[0 0.3 .6]);
set(gca,'tickdir','out',...
	'xlim',[-10e-6 17e-6],...
	'xtick',XV);
title('Winter Curl, Amerasian Basin, 1993-2016');

subplot(2,2,2);
bs=bar(XV,Hs,1);
set(bs,'FaceColor',[0.6 0. 0]);
set(gca,'tickdir','out',...
	'xlim',[-10e-6 17e-6],...
	'xtick',XV);
title('Summer Curl, Amerasian Basin, 1993-2016');

nyr=nc/12;
AA=reshape(CRL,[12,nyr]);
Cm = mean(AA'); % annual mean curl

subplot(2,2,3);
bg=bar(Cm);
set(bg,'FaceColor',[0.5 0.5 0.5]);
set(gca,'tickdir','out',...
	'xlim',[0.5 12.5],...
	'xtick',[1:12]);
title('Vorticity Seasonality');

bottom_text(btx,'pwd',1);

% ---------------
% Plot Fram Wind:
% ---------------
irr=find(TT<1997);
ipp=find(TT>=1997);
mf90=mean(VFr(irr)); % annual mean curl prior 1997
mf00=mean(VFr(ipp)); % annueal mean curl after 1997

iw1=find(TY<1997);
iw2=find(TY>=1997);
wmf90=mean(VFrW(iw1)); % winter mean Fram winds prior 1997
wmf00=mean(VFrW(iw2)); % winter mean Fram winds >1997
smf90=mean(VFrS(iw1)); % winter mean Fram winds prior 1997
smf00=mean(VFrS(iw2)); % winter mean Fram winds >1997


figure(5); clf;
axes('Position',[0.08 0.58 0.85 0.35]);
hold on;
yl1 = -9;
yl2 = 9;
for iyr=1993:2017
  ix1=iyr-0.25;
  ix2=iyr+0.25;
  xx=[ix1 ix2 ix2 ix1];
  yy=[yl1 yl1 yl2 yl2];
  pp = patch(xx,yy,[0.9 0.99 1]);
  set(pp,'Edgecolor','none');
end
plot(TT,VFr,'k');
plot([TT(1) TT(max(irr))],[mf90 mf90],'r--');
plot([TT(min(ipp)) TT(max(ipp))],[mf00 mf00],'g--');

set(gca,'tickdir','out',...
	'xlim',[TT(1) TT(end)],...
	'xtick',[floor(TT(1)):1:ceil(TT(end))],...
	'ylim',[yl1 yl2],...
	'xgrid','on');
title('Wind Fram Strait');



axes('Position',[0.08 0.08 0.85 0.35]);
bb=bar(TY,VFrW);
hold on;
plot([1992.5 1996.25],[wmf90 wmf90],'r--');
plot([1996.7 2016.5],[wmf00 wmf00],'g--');

set(bb,'FaceColor',[0. 0.6 0.8],'EdgeColor','none');
set(gca,'tickdir','out',...
	'xlim',[1992.5 2016.5],...
	'xtick',[1993:2016]);
title('Winter Wind , Fram Strait');

bottom_text(btx,'pwd',1);


figure(6); clf;
axes('Position',[0.08 0.58 0.85 0.35]);
bb=bar(TY,VFrS);
hold on;
plot([1992.5 1996.25],[smf90 smf90],'r--');
plot([1996.7 2016.5],[smf00 smf00],'g--');

set(bb,'FaceColor',[0. 0.6 0.8],'EdgeColor','none');
set(gca,'tickdir','out',...
	'xlim',[1992.5 2016.5],...
	'xtick',[1993:2016]);
title('Summer Wind , Fram Strait');

bottom_text(btx,'pwd',1);


