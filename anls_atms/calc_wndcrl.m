% Calculate monthly wind stress curl over the central Arctic
% and wind over the Fram &  Bering straits
% from NCEP CFSR and CFSv2
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

s_mat = 1;

pthdat1 = '/Net/kronos/ddmitry/ncep_cfsr/';
pthdat2 = '/Net/kronos/ddmitry/ncep_cfsv2/';
pthmat  = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.72/data_mat/';
pth72   = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.72/topo_grid/';

btx = 'calc_wndcrl.m';

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
    
% Mean curl in Beaufort Gyre
    sC=0;
    sA=0;
    ik1=42;
    ik2=120;
    jk1=150;
    jk2=200;
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
clear CW CS CY
for iyr=1993:2016
  cp=cp+1;
  IY = find(DV(:,1)==iyr);
  dmm=mean(CRL(IY));
  CY(cp,1)=dmm; % annual mean
%  IYW=find(DV(:,1)==iyr & (DV(:,2)<4 | DV(:,2)>9));
  IYW=find((DV(:,1)==iyr-1 & DV(:,2)>9) |...
	   (DV(:,1)==iyr & DV(:,2)<=4));
  dmm=mean(CRL(IYW));
  CW(cp,1)=dmm; % winter curl by years
  TY(cp,1)=iyr;
  
  IYS=find((DV(:,1)==iyr & (DV(:,2)>4 & DV(:,2)<10)));
  dmm=mean(CRL(IYS));
  CS(cp,1)=dmm; % summer curl by years
  
end

irr=find(TT<1997);
ipp=find(TT>=1997);
mn90=mean(CRL(irr)); % annual mean curl prior 1997
mn00=mean(CRL(ipp)); % annueal mean curl after 1997

iw1=find(TY<1997);
iw2=find(TY>=1997);
wmn90=mean(CW(iw1)); % winter mean curl prior 1997
wmn00=mean(CW(iw2)); % winter mean curl >1997

WCRL.Title='Wind Curl over Amerasian Basin, CFSR/CFSv2';
WCRL.TM            = DNM;
WCRL.Wind_Curl     = CRL;
WCRL.Wind_Yrs      = TY;
WCRL.Wind_Annual   = CY;
WCRL.Wind_Curl_wnt = CW;
WCRL.Wind_Curl_smr = CS;

if s_mat==1
  fout=sprintf('%swnd_curl_Arctic_cfsr.mat',pthmat);
  fprintf('Saing wind curl %s\n',fout);
  save(fout,'WCRL');
end




%For it=1:nc
figure(3); clf;
axes('Position',[0.08 0.58 0.85 0.35]);
hold on;
yl1 = -9e-6;
yl2 = 5e-6;
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
title('Wind curl, Amerasian Basin');



axes('Position',[0.08 0.08 0.85 0.35]);
bb=bar(TY,CW);
hold on;
plot([1992.5 1996.25],[mn90 mn90],'r--');
plot([1996.7 2016.5],[mn00 mn00],'g--');

set(bb,'FaceColor',[0. 0.6 0.8],'EdgeColor','none');
set(gca,'tickdir','out',...
	'xlim',[1992.5 2016.5],...
	'xtick',[1993:2016]);
title('Winter Wind curl, Amerasian Basin');

bottom_text(btx,'pwd',1);


XV = [-9:1:5]*1e-6;
Hw = hist(CRL(IW),XV);
Hs = hist(CRL(IS),XV);

figure(4); clf;
subplot(2,2,1);
bw=bar(XV,Hw,1);
set(bw,'FaceColor',[0 0.3 .6]);
set(gca,'tickdir','out',...
	'xlim',[-10e-6 6e-6],...
	'xtick',XV);
title('Winter Curl, Amerasian Basin, 1993-2016');

subplot(2,2,2);
bs=bar(XV,Hs,1);
set(bs,'FaceColor',[0.6 0. 0]);
set(gca,'tickdir','out',...
	'xlim',[-10e-6 6e-6],...
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
	'xtick',XV);
title('Vorticity Seasonality');

bottom_text(btx,'pwd',1);

