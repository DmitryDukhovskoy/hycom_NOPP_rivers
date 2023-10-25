% Try assimilate "in situ" observations
% into HYCOM S field 
% within a specifed box 
% - satellite footprint
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;

startup;

close all
clear

regn = 'ARCc0.04';
%expt = 011;  
expt = 012;  

iyr = 2005;

pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';

ftopo = sprintf('%s/depth_%s_17DD.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2


yr=iyr;
iday=180;
%  pthbin = '/Net/mars/ddmitry/hycom/ARCc0.04/output/';
pthbin = sprintf('/nexsan/hycom/ARCc0.04_%3.3i/data/%i/',expt,yr);  
%  pthbin = sprintf('/nexsan/hycom/ARCc0.04_011/data012/%i/',yr);  % Greenland on exp

fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);

dnmb=datenum(yr,1,1)+iday-1;
DV=datevec(dnmb);
mo=DV(2);
mday=DV(3);

plr=1;
[F,n,m,l] = read_hycom(fina,finb,'salin','r_layer',plr);
F(F>1e6)=nan;
SS = squeeze(F);

% Define boxes:
dbx=100e3; % approximate footprint dimensions
BX=[   1360        780
        849        1171
         854        1171
         989         915
        1140         859
        1663        1148
        1659         694
        2017        1536];


nbx=size(BX,1);
for ibb=1:1
  i0=BX(ibb,1);
  j0=BX(ibb,2);
  dx0=DX(j0,i0);
  dy0=DY(j0,i0);
  nx=round(dbx/dx0);
  ny=round(dbx/dy0);
  nx2=round(nx/2);
  ny2=round(ny/2);
  h=HH(j0-ny2:j0+ny2,i0-nx2:i0+nx2);
  inn=find(isnan(h));
  if ~isempty(inn); 
    fprintf('ibb=%i, land points in the box Nland=',ibb,length(inn));
  end
  
  S0=SS(j0-ny2:j0+ny2,i0-nx2:i0+nx2);
  [ms,ns]=size(S0);
  ntot=ms*ns;

%
% Generate N observations:
  Nobs=round(0.1*ntot);
%  IIo=round(rand(Nobs,1)*(ns-1)+1);
% Random point in the upper right quadrant
  IIo=round(rand(Nobs,1)*(ns/2-1)+ns/2);
  JJo=round(rand(Nobs,1)*(ms/2-1)+ms/2);

  [ha,hx]=hist(S0,30);
  dhx=hx(2)-hx(1);
%  han=ha/(dhx*sum(ha)); % probability
  han=ha/length(S0);  % frequency

% Kalman filtering  
% X(k) = A*X(k-1) + w(k-1) - process
  Xapr=S0(:); % x apriori
  nk=length(Xapr);
  
% Create Matrix A:
  A=eye(nk);
  
% Process noise covariance
% assume constance covariance - this 
% will be estimated from HYCOM
% X(k) = A*X(k-1) + w(k-1)
% Q=E[w(k-1)*w'(k-1)]
  sgm2=0.1;       % variance 
  Q=sgm2*eye(nk);   % process noise covariance
  
% Observations: replace model data with obs
% that need to assimilate
% Zk = H*X(k)+nu(k)
  Iobs=sub2ind(size(S0),JJo,IIo);       % indices for observations
  mk=length(Iobs);
  Sobs=Xapr(Iobs)+rand*0.3;  % observed values
  Zk=Sobs(:);           % mx1 matrix
  R=0.0001*eye(mk);  % measurement noise covariance
  
% Matrix H (m x n): relates model to obs
% assumes perfect relation with no noise
  H=zeros([mk,nk]);
  for ib=1:mk
    ii0=Iobs(ib);
    H(ib,ii0)=1;
  end

% Initial esimates:  
% A priori estimate error covariance
% start with something
  Eapr=Xapr*0+0.01;
  Ppst=Eapr*Eapr';   % P(k) apriori=E[e(k)*e(k)'] initial P(k-1)
  Xpst=Xapr;
  
  for ikk=1:2  % time update - several iterations
% Time Update - predict:    
% (1) Project the state ahead:
% Xapr(k)=A*Xpst(k-1)+ ... - in our case - steady state, no projection
   Xapr=Xpst;
% (2) Project the error covariance ahead:
% P(k) = A*P(k-1)*A'+Q
    Papr = A*Ppst*A'+Q;
    
% Measurement Update (Correct):
% Step (1) update Gain:  
% Kalman Gain:
    K = Papr*H'/(H*Papr*H'+R);
%
% (2) Update an estimate with observation:
    Xpst = Xapr + K*(Zk - H*Xapr); 
    
% (3) Update Covariance:
    Ppst = Papr - K*H*Papr;
  end
  
end

btx='spatialS_KalmanF.m';

figure(3); clf;
plot(S0(:));
hold
plot(Iobs,Sobs,'r*');
plot(Xpst,'g');


c1=34.5;
c2=34.9;
CMP=create_colormap7(400,c1,c2);
cmp=CMP.colormap;


figure(1); clf;
pcolor(S0); shading flat;
hold on;
plot(IIo,JJo,'k.','Markersize',14);
caxis([c1 c2]);
colormap(cmp);
hb=colorbar;
stl=sprintf('S from HYCOM, i=%i, j=%i',i0,j0);
title(stl);

%btx='spatialS_KalmanF.m';
bottom_text(btx,'pwd',1);


Sflt=reshape(Xpst,[ms,ns]);
figure(2); clf;
pcolor(Sflt); shading flat;
caxis([c1 c2]);
colormap(cmp);
hb=colorbar;
title('Assimilated S');

btx='spatialS_KalmanF.m';
bottom_text(btx,'pwd',1);


% Histograms
ss1=S0(:);
%hist(S0(:));
[ha,hx]=hist(ss1(:),30);
dhx=hx(2)-hx(1);
han=ha/length(ss1);  % frequency
Smn=mean(ss1);
sg1=Smn-std(ss1);
sg2=Smn+std(ss1);

figure(4); clf;
axes('Position',[0.1 0.6 0.8 0.3]);
hbb=bar(hx,han);
hold on;
plot([Smn Smn],[0 max(han)],'r--');
plot([sg1 sg1],[0 max(han)],'k--');
plot([sg2 sg2],[0 max(han)],'k--');

% Plot mean estimates from random "observations":
Smobs=mean(S0(Iobs));
plot([Smobs Smobs],[0 max(han)],'b--');


xl1=min(hx)-dhx/2;
%xl2=max([max(Xpst),max(S0(:))]);
xl2=max(ss1)+dhx/2;
set(hbb,'Facecolor',[0.6 0.6 0.6]);
set(gca,'tickdir','out',...
	'xtick',hx,...
	'xlim',[xl1 xl2]);

stl=sprintf('i=%i, j=%i, S, hycom04-012, %s',i0,j0,datestr(dnmb));
title(stl);
bottom_text(btx,'pwd',1,'Position',[0.05 0.4 0.4 0.05]);


% Assimilated data:
ssf=Sflt(:);
[hfa,hfx]=hist(ssf,30);
dfhx=hfx(2)-hfx(1);
hfan=hfa/length(Sflt);  % frequency
Sfmn=mean(ssf);
sfg1=Sfmn-std(ssf);
sfg2=Sfmn+std(ssf);

axes('Position',[0.1 0.15 0.8 0.3]);
hbb=bar(hfx,hfan);
hold on;
plot([Sfmn Sfmn],[0 max(hfan)],'r--');
plot([sfg1 sfg1],[0 max(hfan)],'k--');
plot([sfg2 sfg2],[0 max(hfan)],'k--');

set(hbb,'Facecolor',[0.6 0.6 0.6]);
set(gca,'tickdir','out',...
	'xtick',hfx,...
	'xlim',[xl1 xl2]);

stl=sprintf('KalmFltr, i=%i, j=%i, S',i0,j0);
title(stl);

%btx='spatialS_KalmanF.m';
bottom_text(btx,'pwd',1);


  