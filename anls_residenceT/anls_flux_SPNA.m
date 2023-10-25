% To test hypothesis about convergence/divergence
% of particles in the SPNA at different depths
%
% Estimate divergence - volume averaged
% using Gauss (diverg) theorem: integrate fluxes over the 
% boundaries
% Time series extracted in calc_flux_SPNA.m
%
addpath /Net/Movies0/ddmitry/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /home/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers;
%addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
%addpath /Net/Movies0/ddmitry/MyMatlab/seawater
startup

close all
clear


pthmat = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/110/data_GG_prt/';

btx = 'anls_flux_SPNA.m'; 

% Mean fluxes by layers
YR1=2003;
YR2=2006;
tcc=0;
for YR=YR1:YR2
  fmat = sprintf('%shycom008_112_SPNAcntr_VolFlux_%i.mat',pthmat,YR);
  fprintf('Loading %s\n',fmat);
  load(fmat);

% Integrate by layers over the whole contour surface
  tmm = GVHFLX.Time;
  nt=length(tmm);
%
% Calculate means and adjust small biases to have overall transport = 0
  cc=0;
  for it=1:nt
    VV  = squeeze(GVHFLX.VflxPnts(it,:,:));
    sv = nansum(VV,2);

    dzz = squeeze(GVHFLX.DZ(it,:,:));
    dzz(dzz==0)=nan;
%    dV=VV./dzz; % flux m2/s or m3/s per 1 m depth
%    dV(isnan(dV))=0;
%    sumV = sum(dV,2);
    dzm = nanmean(dzz,2);
%    zzm = -cumsum(dzm);

    cc=cc+1;
    tcc=tcc+1;
    MNV.TM(tcc)=tmm(cc);
%    MNV.VFlx_m2s(:,tcc)=sumV; % integrated around contour Vol flx by layers / dz = m2/s
%    MNV.Zm(:,tcc)=zzm;
    MNV.VFlx_m3s(:,tcc)=sv;  % integrated over the lateral surfaces
    MNV.dZm(:,tcc)=dzm;
  end

end


%vfl = MNV.VFlx_m2s;
%mnV = mean(vfl,2);

% Interpolate onto fine res vert grid
dzi = 2;
ZZi = [0:-dzi:-4000];

% Winter/summer mean divergence vol integrated
% by layers m3/s
TM = MNV.TM;
DV = datevec(TM);

Is = find(DV(:,2)>=5 & DV(:,2)<=9);
Iw = find(DV(:,2)>=11 | DV(:,2)<=3);

% Winter
Vw  = MNV.VFlx_m3s(:,Iw);
dzw = MNV.dZm(:,Iw);
nVw = Vw./dzw; % m2/s
mVw = nanmean(nVw,2); % mean transp m2/s
mVw(end) = mVw(end-1); % get rid of nans for interpolation
mVw=[mVw(1);mVw];      % for interpolation
mDZ = nanmean(dzw,2);
mDZ(41)=10;  % for interpolation extend bottom
Zw = -cumsum(mDZ);
Zw = [0;Zw];

mVwi = interp1(Zw,mVw,ZZi,'pchip');
mVwi = - mVwi; % make div>0 = divergence
vtot = sum(mVwi);
% To balance 
dv=-vtot/length(mVwi);
mVwi=mVwi+dv;
dvw=dv;

[a1,a2]=size(nVw);

clear pw1 pw9
for k=1:a1
  pw1(k,1)=prctile(nVw(k,:),10);
  pw9(k,1)=prctile(nVw(k,:),90);
end 

pw1 = [pw1(1);pw1];
pw1(end) = pw1(end-1);
pw1i = interp1(Zw,pw1,ZZi,'pchip');
 
pw9 = [pw9(1);pw9];
pw9(end) = pw9(end-1);
pw9i = interp1(Zw,pw9,ZZi,'pchip');

pw1i=pw1i+dv;
pw9i=pw9i+dv;



% Summer
Vs  = MNV.VFlx_m3s(:,Is);
dzs = MNV.dZm(:,Is);
nVs = Vs./dzs; % m2/s
mVs = nanmean(nVs,2); % mean transp m2/s
mVs(end) = mVs(end-1); % get rid of nans for interpolation
mVs=[mVs(1);mVs];      % for interpolation
mDZ = nanmean(dzs,2);
mDZ(41)=10;  % for interpolation extend bottom
Zs = -cumsum(mDZ);
Zs = [0;Zs];

mVsi = interp1(Zs,mVs,ZZi,'pchip');
mVsi = -mVsi; % make div>0 = divergence
vtot = sum(mVsi);
% To balance 
dv=-vtot/length(mVsi);
mVsi=mVsi+dv;
dvs=dv;

% Overall mean:
Vy  = MNV.VFlx_m3s;
dzy = MNV.dZm;
nVy = Vy./dzy; % m2/s
mVy = nanmean(nVy,2); % mean transp m2/s
mVy(end) = mVy(end-1); % get rid of nans for interpolation
mVy=[mVy(1);mVy];      % for interpolation
mDZ = nanmean(dzs,2);
mDZ(41)=10;  % for interpolation extend bottom
Zy = -cumsum(mDZ);
Zy = [0;Zy];

mVyi = interp1(Zy,mVy,ZZi,'pchip');
mVyi = -mVyi; % make div>0 = divergence
vtot = sum(mVyi);
% To balance 
dv=-vtot/length(mVyi);
mVyi=mVyi+dv;


% Plot histograms of div at depth levels where Lagr particles 
Ik=[10;15;23;31];
nL=length(Ik);
for jk=1:nL
  iz=Ik(jk);
  dmm = -(nVw(iz,:)+dvw); % invert sign --> div>0 - divergence
  [N,X]=hist(dmm,20);
  HSTW(jk).N=N;
  HSTW(jk).X=X;

  dmm = - (nVs(iz,:)+dvs);
  [N,X]=hist(dmm,20);
  HSTS(jk).N=N;
  HSTS(jk).X=X;
end

    
% Plot winter/summer vol integrated divergences by layers
% Units are m3/s per 1 m layer thickness (i.e. m2/s)
figure(1); clf;
axes('position',[0.09 0.1 0.4 0.8]);
hold on;
plot(mVyi,ZZi,'linewidth',2.5,'Color',[0 0 0]); 
plot(mVwi,ZZi,'linewidth',2.5,'Color',[0 0.4 0.8]);  
plot(mVsi,ZZi,'linewidth',2.5,'Color',[0.8 0.3 0]);  

%plot(-pw1i,ZZi);
%plot(-pw9i,ZZi);
set(gca,'tickdir','out',...
        'xlim',[-5e4 10e4],...
        'ylim',[-4000 0],...
        'xtick',[-6e4:2e4:10e4],...
        'ytick',[-4000:250:0],...
        'xgrid','on',...
        'ygrid','on',...
        'fontsize',14);
lgd = legend('Overall','Winter','Summer');
set(lgd,'Position',[0.65 0.80 0.15 0.12],'Fontsize',14);
stl = sprintf('VolIntgr Div(U)/1m (m2/s), div>0/conv<0, ARCc008, %i-%i',YR1,YR2);
title(stl);

% Zoom in upper layers
axes('Position',[0.55 0.1 0.4 0.4]);
hold on;
plot(mVyi,ZZi,'linewidth',2.5,'Color',[0 0 0]);
plot(mVwi,ZZi,'linewidth',2.5,'Color',[0 0.4 0.8]);
plot(mVsi,ZZi,'linewidth',2.5,'Color',[0.8 0.3 0]);

%plot(-pw1i,ZZi);
%plot(-pw9i,ZZi);
set(gca,'tickdir','out',...
        'xlim',[-5e4 10e4],...
        'ylim',[-600 0],...
        'xtick',[-6e4:2e4:10e4],...
        'ytick',[-1000:50:0],...
        'xgrid','on',...
        'ygrid','on',...
        'fontsize',14);

bottom_text(btx,'pwd',1);


% Cumulative transport acros the lateral surfaces 
cVyi=cumsum(mVyi)*dzi*1e-6;  % Sv
cVwi=cumsum(mVwi)*dzi*1e-6;
cVsi=cumsum(mVsi)*dzi*1e-6;

figure(2); clf;
axes('position',[0.09 0.1 0.4 0.8]);
hold on;
plot(cVyi,ZZi,'linewidth',2.5,'Color',[0 0 0]);
plot(cVwi,ZZi,'linewidth',2.5,'Color',[0 0.4 0.8]);
plot(cVsi,ZZi,'linewidth',2.5,'Color',[0.8 0.3 0]); 

%plot(-pw1i,ZZi);
%plot(-pw9i,ZZi);
set(gca,'tickdir','out',...
        'xlim',[-6 6],...
        'ylim',[-4000 0],...
        'xtick',[-6:6],...
        'ytick',[-4000:250:0],...
        'xgrid','on',...
        'ygrid','on',...
        'fontsize',14);
lgd = legend('Overall','Winter','Summer');
set(lgd,'Position',[0.65 0.80 0.1 0.08],'Fontsize',12);
stl = sprintf('Depth-intgr, Net Vol flux across lat bound (+ out), Sv, ARCc008, %i-%i',YR1,YR2);
title(stl);

bottom_text(btx,'pwd',1);







