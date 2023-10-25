% Calculate & plot 
% heat flux through major straits 
% See Schauer et al., 2005, JGR Arctic warmin gthrough the Fram Strait
% Uses Tref = -0.1C (Similar to Mosby 1962; Aagaard and Greisman, 1975; and -0.7C)
% Better use approach of "pipe" - discussed in later paper Schauer 
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

Cp = 4.184*1e3;  %J/kg*C water heat capacity
rho = 1025;      % kg/m3
Tref = -0.1; % Reference T for calculating heat flux

expt = '110';
TV   = 11;  % topo version
SEGM{1} = 'BeringS';
SEGM{2} = 'FramS';
SEGM{3} = 'BarOp';
SEGM{4} = 'DavisS';
nsgm = length(SEGM);
txb = 'plot_heat_flx.m';

s_mat=0;
%Zmn = 100; % average over the top Zmn m
rg  = 9806;
hgg = 1e20; % 

pthmat  = '/Net/mars/ddmitry/hycom/ARCc0.08/data_mat/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
%monmat = sprintf('%sEKE_meanUV_lev%i.mat',pthmat,Nlev);

% T/S sections in straits
fmat = sprintf('%s%s_TS_straits_%4.4i_%4.4i.mat',...
	       pthmat,expt,1993,2016);
fprintf('Loading %s\n',fmat);
load(fmat); % SGMTS

% UV sections in straits
fmat = sprintf('%s%s_UV_straits_%4.4i_%4.4i.mat',...
	       pthmat,expt,1993,2016);
fprintf('Loading %s\n',fmat);
load(fmat);  % UV

% Time - time series can be of different lenghts
TMts= SGMTS(1).TM;
nTS = length(TMts); 
d1=TMts(1);
dv=datevec(d1);
DVts = datevec(TMts);
yr1=dv(1);
dj1=datenum(yr1,1,1);

TMuv= UV(1).TM;
nUV = length(TMuv);
DVuv= datevec(TMuv);

% Match time series:
cc=0;
for ik=1:min([nTS,nUV]);
  d1=TMts(ik);
  i2 = find(TMuv==d1);
  if isempty(i2), continue; end;
  cc = cc+1;
  Tindx(cc,1) = ik;
  Tindx(cc,2) = i2;
end


% Start from:
ftopo = sprintf('%sdepth_ARCc0.08_%2.2i.nc',pthtopo,TV); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');

[mm,nn]=size(LON);

for isgm=1:nsgm
  segm=SEGM{isgm};
  fprintf('Calculating flux for %s\n',segm);
  
  switch(segm)
   case('BeringS');
    is1=634;
    js1=1919;
    is2=658;
    js2=js1;
    X=LON(js1,is1:is2);
    Y=LAT(js1,is1:is2);
    I=[is1:is2]';
    J=ones(size(I))*js1;
   case('FramS');
    is1=935;
    js1=959;
    is2=1070;
    js2=js1;
    X=LON(js1,is1:is2);
    Y=LAT(js1,is1:is2);
    I=[is1:is2]';
    J=ones(size(I))*js1;
   case('BarOp');
    is1=1166;
    js1=456;
    is2=is1;
    js2=929;
    X=LON(js1:js2,is1);
    Y=LAT(js1:js2,is1);
    J=[js1:js2]';
    I=ones(size(J))*is1;
   case('DavisS');
    is1=478;
    js1=680;
    is2=568;
    js2=680;
    X=LON(js1,is1:is2);
    Y=LAT(js1,is1:is2);
    I=[is1:is2]';
    J=ones(size(I))*js1;
  end

  p_map=0;
  if p_map>0
    figure(11); clf;
    contour(HH,[0 0],'k');
    hold on
    plot(I,J,'r.-');
%    keyboard
  end


  dmm = SGMTS(isgm).Temp;
  Ft  = dmm(Tindx(:,1),:,:);
  dmm = UV(isgm).UV_normal;
  Fuv = dmm(Tindx(:,2),:,:);
  dmm = UV(isgm).ZZ;
  ZZ  = dmm(Tindx(:,2),:,:);
  nrc = size(Fuv,1);
  dX  = UV(isgm).Dist;
  
  TFLX(isgm).segm = segm;
  TFLX(isgm).TM = TMuv(Tindx(:,1));
  
  for jn=1:nrc
    zz = squeeze(ZZ(jn,:,:));
    dZ = abs(diff(zz,1));
    ll  = size(dZ,1);
    uv = squeeze(Fuv(jn,:,:));
    tt = squeeze(Ft(jn,:,:));
% Limit only for t>Tref
    Irf = find(tt<Tref);
    Iatl = find(tt<1); % Atlantic water only > 1C
    Ipls = find(uv>=0); % positive flow (can be in or out)
    Imns = find(uv<0);
    
    Acell = dZ*0;
    for k=1:ll
      Acell(k,:) = dZ(k,:).*dX';
    end
    
    Flx = rho*Cp*(tt-Tref).*uv; % Flux J/s*m2
    tFlx = nansum(nansum(Flx.*Acell)); % J/s through strait
    tFlxPl= nansum(nansum(Flx(Ipls).*Acell(Ipls))); % flux with U>0
    tFlxMn= nansum(nansum(Flx(Imns).*Acell(Imns))); % U<0
    
    
    dmm = Flx;
    dmm(Irf)=nan; % only for T>-0.1
    tFlxR = nansum(nansum(dmm.*Acell)); % J/s through strait
    dmm = Flx;
    dmm(Iatl)=nan; % only for T>-0.1
    tFlxA = nansum(nansum(dmm.*Acell)); % J/s through strait
    
    f_plt=0;
    if f_plt==1
      XX = cumsum(dX)-dX(1);
      dmm = Flx;
      dmm = [Flx(1,:);dmm];
      dmm = dmm*1e-12;  % TJ
      
      figure(1); clf;
      pcolor(XX,zz,dmm); shading flat;
      title(sprintf('Heat Flux, TJ/s, %s',segm));
    end;
    
    TFLX(isgm).TFlx_dpth(jn)= tFlx;  % total depth with no other constraints
    TFLX(isgm).TFlx_trf(jn) = tFlxR; % only for T>Tref
    TFLX(isgm).TFlx_atl(jn) = tFlxA; % only for T>1C
    TFLX(isgm).TFlx_Pl(jn)  = tFlxPl; % Tflx with U>0
    TFLX(isgm).TFlx_Mn(jn)  = tFlxMn; % Tflx with U<0
    
  end
end

  
isgm=2;
segm = TFLX(isgm).segm;
TM  = TFLX(isgm).TM;
DV  = datevec(TM(1));
TMd = TM-TM(1); % days;
Tyr=TMd/365.25+DV(1,1);

cf=1e-12;
ff = TFLX(isgm).TFlx_dpth*cf; % TJ/
ffr = TFLX(isgm).TFlx_trf*cf; % TJ/
ffp = TFLX(isgm).TFlx_Pl*cf; 
ffm = TFLX(isgm).TFlx_Mn*cf;
% Butterworth filter:
% !!! Note in Matlab: cutoff freq. has to be normalized
% so that it is 0<fcutoff<1,
% with 1 corresponding to half the sample rate !!!
Dflt = 180;    % flt freq. 
dday=5;        % output freq.
sfr =1/(dday); % sampl. freq: 5 days
fcutoff=1/Dflt; 
fnorm = fcutoff/(sfr/2);
[Bt1,Bt2]=butter(5,fnorm,'low');  % low-pass filter
Flow=filtfilt(Bt1,Bt2,ff);  % frwd/reverse filtering
ffl=Flow;
Flow=filtfilt(Bt1,Bt2,ffr);  % frwd/reverse filtering
ffrl=Flow;
Flow=filtfilt(Bt1,Bt2,ffp);  % frwd/reverse filtering
ffpl=Flow;
Flow=filtfilt(Bt1,Bt2,ffm);  % frwd/reverse filtering
ffml=Flow;


figure(2); clf;
axes('position',[0.08 0.6 0.85 0.32]);
plot(Tyr,ff,'Color',[0.8 0.8 0.8]);
hold on;
plot(Tyr,ffl,'r');
set(gca,'xlim',[floor(Tyr(1)) ceil(Tyr(end))],...
	'tickdir','out',...
	'xtick',[round(Tyr(1)):ceil(Tyr(end))+1]);

mn=mean(ffl);
sgm=std(ffl);
stl=sprintf('NetTot HFlx, Tref=%3.2fC, <F>=%4.2f+/-%4.2f TJ, %s',Tref,mn,sgm,segm);
title(stl);

axes('position',[0.08 0.08 0.85 0.32]);
plot(Tyr,ffr,'Color',[0.8 0.8 0.8]);
hold on;
plot(Tyr,ffrl,'r');
set(gca,'xlim',[floor(Tyr(1)) ceil(Tyr(end))],...
	'tickdir','out',...
	'xtick',[round(Tyr(1)):ceil(Tyr(end))+1]);

mn=mean(ffr);
sgm=std(ffr);
stl=sprintf('Net HFlx (T>Tref), <F>=%4.2f+/-%4.2f TJ, %s',mn,sgm,segm);
title(stl);
bottom_text(txb,'pwd',1);

% Plot in/out flows
figure(3); clf;
axes('position',[0.08 0.6 0.85 0.32]);
plot(Tyr,ffp,'Color',[0.8 0.8 0.8]);
hold on;
plot(Tyr,ffpl,'r');
set(gca,'xlim',[floor(Tyr(1)) ceil(Tyr(end))],...
	'tickdir','out',...
	'xtick',[round(Tyr(1)):ceil(Tyr(end))+1]);

mn=mean(ffp);
sgm=std(ffp);
stl=sprintf('HFlux U>0, <F>=%4.2f+/-%4.2f TJ, %s',mn,sgm,segm);
title(stl);

axes('position',[0.08 0.08 0.85 0.32]);
plot(Tyr,ffm,'Color',[0.8 0.8 0.8]);
hold on;
plot(Tyr,ffml,'r');
set(gca,'xlim',[floor(Tyr(1)) ceil(Tyr(end))],...
	'tickdir','out',...
	'xtick',[round(Tyr(1)):ceil(Tyr(end))+1]);

mn=mean(ffm);
sgm=std(ffm);
stl=sprintf('HFlux U<0, <F>=%4.2f+/-%4.2f TJ, %s',mn,sgm,segm);
title(stl);
bottom_text(txb,'pwd',1);



