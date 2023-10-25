% Analyze GFWA fluxes in the SPNA
% majors straits and SPNA contour
% data extracted in extr_Trcr_daily_straits08.m
% Tracer is converted into GFWA
% using fract of Tr in a cell given overall Tracer mass in the domain

addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/Wavelet;
startup;

close all
clear


expt=112;
TV=11;
YR1=1995;
YR2=2016;
dday=7;
nTr=1;  % extract 1 tracer at a time, 1 - Greenland

s_mat=0;  % ==2 - start from last saved rec # in YEAR

hgg=1e20;
f_zgrd=0;  % =1 - calculate fluxes from z-grid interpolated U,T,S - less accurate
           % mostly for comparison and validation
ptharch = '/nexsan/people/ddmitry/';
pthout  = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_straits/';
pthmat  = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_theresa/';
pthmat3 = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
pthriv = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/110/data_mat/';
pthtopo= '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';

btx='anls_TrFlux_SPNA.m';

fprintf('arc08-%3.3i Tracer Fluxes Gates %i-%i, save=%i\n',...
        expt,YR1,YR2,s_mat);


ftopo = sprintf('%sdepth_ARCc0.08_%2.2i.nc',pthtopo,TV); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2


% Fram Section is close to Moorings ~79N
SCT = sub_define_SPNA_sections(HH,LON,LAT);
nsct = length(SCT);


CLR=[0 0.4 0.8; ...
     0.8 0.4 0; ...
     1 0.2 0; ...
     0 1 0; ...
     0.8 0 0.6; ...
     0.8 1 0; ...
     0.3 0.7 1; ...
     1 0.7 0.5; ...
     0.9 0.3 0.6; ...
     0.2 0.7 0.4; ...
     0.1 0.4 0.3; ...
     0.8 0.3 0.9; ...
     1 0.4 0.6];

f_map=0;
if f_map==1
  fprintf('Drawing map with segments\n');
  fn=1;
%  sub_plot_Greenl_contour(HH,LON,LAT,fn,GC);
  figure(1); clf;
  contour(HH,[0 0],'k');
  hold on;
  contour(HH,[-5000:500:-100],'Color',[0.9 0.9 0.9]);

  for ip=1:nsct
%    IJ=SCT(ip).IJ;
%    plot([IJ(1,1) IJ(2,1)],[IJ(1,2) IJ(2,2)],...
%        'Linewidth',2.5,'Color',[1. 0.6 0]);
   clr=CLR(ip,:);
    IIs=SCT(ip).I;
    JJs=SCT(ip).J;
    plot(IIs,JJs,'-',...
         'Linewidth',2.5,'Color',clr);
    Ip=SCT(ip).IJPR(1);
    Jp=SCT(ip).IJPR(2);

    plot(Ip,Jp,'.','Markersize',14,'Color',clr);
  end

  axis('equal');
  set(gca,'xlim',[300 1300],...
          'ylim',[100 1300]);

  bottom_text(btx,'pwd',1);

%keyboard
end

frv = sprintf('%sGreenland_cumFWFlux.mat',pthriv);

f_griv = 0; %=1 - rederive cumulative Greenland FWFlux from Bamber
if f_griv==1
  sub_get_GrRunoff(friv)
end
fprintf('f_griv %i, Loading %s\n',f_griv,frv);
load(frv);


% Get tracer fraction in grid cells
IntSrf = 1; % = 1 - integrates over all layers from 0 - ilv, for ilv=1 is the same
%ilv = 1; % 0-50m
%ilv = 2; % 50-150m
%ilv = 3; % 150-300m
%ilv = 4; % 300-500 m
ilv = 5; % whole depth <- do not use this

% Vertical layers
LRS = load('LRS.dat');
nlrs= length(LRS);


zz1 = LRS(ilv,1);
zz2 = LRS(ilv,2);

if IntSrf==1; zz1=0; end;

dz=abs(zz2-zz1);
hZ = abs(LRS(ilv,2)-LRS(ilv,1));

fprintf('FW Volume Subpolar Gyre, ilv=%i, %i - %i, nTr=%i, 1993-2016\n',ilv,zz1,zz2,nTr);


cc=0;
for YR=YR1:YR2
  fmatout=sprintf('%shycom008_%3.3i_Trcr%2.2i_StraitFluxesDay_%4.4i.mat',...
                    pthmat,expt,nTr,YR);
  fprintf('Loading %s \n',fmatout);
  load(fmatout);

  cc=cc+1;
  TM=SCT(1).Time;
  DV=datevec(TM);
  Is=find(DV(:,2)>4 & DV(:,2)<10);
  Iw=find(DV(:,2)<=4 | DV(:,2)>=10);

%
% Obtain time series of overall Mass tracer for all time in this year:
  TM = SCT(1).Time;
  DVm = datevec(TM);
  imo0=-1;
  MTr_dom=[];
  for it=1:length(TM)
    DV = DVm(it,:);
    iyr = DV(1);
    imo = DV(2);
    if imo~=imo0
      expt0=110;
      %pthmat2=sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt0);
      pthmat2=sprintf('/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/%3.3i/data_matTr/',expt);
      fmat2 = sprintf('%sMassTr%2.2i_lrs_%i%2.2i.mat',pthmat2,nTr,iyr,imo);
      fprintf('Loading %s\n',fmat2);
      load(fmat2);
      imo0=imo;
    end
% find whole-depth layer
    nlr = length(TRCR);
    ibtm=5; % whole-depth tracer mass
    Tr_dom=squeeze(TRCR(ibtm).MassTr_kg); % tracer integrated ove whole water depth
    amm=Tr_dom./(abs(HH).*DX.*DY);  % kg/m3
    Tr_dom(amm<=0.125)=nan;
    MTr_dom(it,1) = nansum(nansum(Tr_dom)); % overall mass of the tracer in the domain
  end

%
% Get tracer fluxes across straits
  nsct = length(SCT);
  for isc=1:nsct
    nm=SCT(isc).Name;
    Nrm=SCT(isc).Norm;
    Xn=Nrm(:,1);
    Yn=Nrm(:,2);
%     
    TFlx=SCT(isc).TrFlx;  % kg/s of tracer
%
% Convert kg/s of tracer into m3/s of GFWA:
% Greenland runoff anomaly for given date
    ism = find(Ygr==DV(1,1));
    if isempty(ism), ism=length(cFWF); end; % after 2016 - same 
    fwf0 = cFWF(ism); % km3

    MTr_dom(MTr_dom==0)=nan;
    rr=TFlx./MTr_dom;
    vflx=rr*fwf0;      % km3/sec
    nrc=length(vflx);
    FLX(isc).Name = nm;
    FLX(isc).GrFlx_mean(cc,1)=sum(vflx)/nrc*3600*24*365; % km3/yr
%
% Compute mean vol flux, Sv
    ZZ=SCT(isc).ZZintrp;
    nlr=length(ZZ);
    dZ=abs(diff(ZZ));
    dZ(nlr)=dZ(nlr-1);
    dl=SCT(isc).segm_dL;
    [DL,DZ]=meshgrid(dl,dZ);

    vflx = SCT(isc).VolFlx_m3s;
    VFlx = mean(vflx,1);
    FLX(isc).VolFlux_Sv(cc,1) = VFlx*1e-6; % Sv 
% 
% Negative and positive Tracer fluxes separately:
% Fluxes over whole sections
    [a1,a2,a3] = size(SCT(isc).Unrm);
    vpos = [];
    vneg = [];
    for iday = 1:a1
      un   = squeeze(SCT(isc).Unrm(iday,:,:));
      tn   = squeeze(SCT(isc).Trcr(iday,:,:)); % kg/m3
      vf0  = un.*tn.*DL.*DZ;
      Ipos = find(vf0>0);
      Ineg = find(vf0<0);
      vpos(iday,1) = nansum(vf0(Ipos));
      vneg(iday,1) = nansum(vf0(Ineg));
    end;
% Convert Tracer conc flux to km3/s of GFWA 
    rr = vpos./MTr_dom;
    TFlx_pos = rr*fwf0;  % km3/sec

    rr = vneg./MTr_dom;
    TFlx_neg = rr*fwf0;  % km3/sec

    FLX(isc).GrFlx_neg(cc,1)=sum(TFlx_neg)/nrc*3600*24*365;
    FLX(isc).GrFlx_pos(cc,1)=sum(TFlx_pos)/nrc*3600*24*365;
    
%
  end

end

% Report
fprintf('GFWA fluxes through SPNA straits\n');
fprintf('Positive flux is north and east\n');
for isc=1:nsct
  nm=SCT(isc).Name;
  flx=FLX(isc).GrFlx_mean;
  mflx=mean(flx);
  nn=length(flx);
  serr=std(flx)/sqrt(nn);

  vflx=FLX(isc).VolFlux_Sv;
  mVflx=mean(vflx);
  sVerr=std(vflx)/sqrt(nn);

  fprintf('%2i %s: GFWA    %6.4g +/- %6.4g km3/yr\n', isc,nm,mflx,serr);
  fprintf('        VolFlx  %6.4g +/- %6.4g Sv\n',mVflx,sVerr);
end

%
% Plot time series of the net flux and out-flux across the bndry:
Fsm=[];  % net flux that is <0
Fneg=[]; % pure negative flux - outflow
Fpos=[];
Fout=0; % outflow of the GFWA
for ii=1:nsct
  if ii==3 | ii==4 | ii==5; continue; end;
  sgn=1;
%  if ii==10; sgn=0; end;

  ff = sgn*FLX(ii).GrFlx_mean;
  ngg= FLX(ii).GrFlx_neg;  
  pss= FLX(ii).GrFlx_pos;

  if isempty(Fsm);
    Fsm=ff;
    Fneg=ngg;
    Fpos=pss;
  else
    Fsm=Fsm+ff;
    Fneg=Fneg+ngg;
    Fpos=Fpos+pss;
  end

  if ff<0
    Fout=Fout+ff;
  end 
end
%Fout(end-1)=-190;
Fout(end)=-184; % impact of southern Open BC - tracer is flowing back

Tyr = [YR1:YR2];

figure(1); clf;
axes('Position',[0.09 0.5 0.85 0.42]);
plot(Tyr,Fout);
title(sprintf('expt=%i, GFWA net outflow, km3/yr',expt));
btx = 'anls_TrFlux_SPNA.m'; 
bottom_text(btx,'pwd',1,'position',[0.05 0.4 0.4 0.04]);

f_out=0;
if f_out==1
		figure(2); clf;
		axes('Position',[0.09 0.5 0.85 0.42]);
		plot(Tyr,Fneg+Fpos);
		title(sprintf('expt=%i, GFWA outflow, km3/yr',expt));
		btx = 'anls_TrFlux_SPNA.m';
		bottom_text(btx,'pwd',1,'position',[0.05 0.4 0.4 0.04]);

		figure(3); clf;
		axes('Position',[0.09 0.5 0.85 0.42]);
		plot(Tyr,Fpos);
		title(sprintf('expt=%i, GFWA inflow, km3/yr',expt));
		btx = 'anls_TrFlux_SPNA.m';
		bottom_text(btx,'pwd',1,'position',[0.05 0.4 0.4 0.04]);
end

%
% Gr FW in the domain
% prepared in dS_FWC_timeseries_SubpolarGyre.m
ilv=5;
fmap='SPG_noNorth'; % SPG no North Sea and no Baffin
fextr=sprintf('%sFWC_TimeSer_%s_ilv%i.mat',pthmat3,fmap,ilv);
load(fextr);
[m,n]=size(VFW);
VFWy=mean(VFW)';  % 
VFWy=[10;20;VFWy];  % missing first 2 years from 110 expt
%VF=reshape(VFW,[m*n,1]);
%Tyr2=[YR1:2016];
Tyr2=[1993:2016];

% Estimate inv. residence time scale:
nFout = -Fout; % swap sign of the outflow for comparison
nFout = [0;1;nFout];
kk = nFout./VFWy;

% Same from regression:
X=VFWy*0+1;
X=[X,VFWy];
Y=nFout;
[B,bint,R,rint,stats] = regress(Y,X);
Yft = X*B;

% Plot FWC and outflow:
figure(2); clf;
axes('Position',[0.09 0.5 0.82 0.42]);
hold on;
plot(Tyr2,nFout,'Linewidth',2); % net outflow
plot(Tyr2,Yft,'Linewidth',2,'Color',[0.6 0.85 0.9]);
yl2=200;
set(gca,'tickdir','out',...
        'xtick',[1993:2016],...
        'ytick',[0:50:200],...
        'ylim',[0 200],...
        'xlim',[1993 2016],...
        'Fontsize',14,...
        'xgrid','on',...
        'ygrid','on');
stl = sprintf('%i, GFWA in SPNA (m3) and net GFWA outflow km3/yr',expt);
title(stl);

axes('Position',[0.09 0.5 0.82 0.42]);
plot(Tyr2,VFWy,'Color',[0.7 0.3 0],'Linewidth',2); % GFWA in SPNA
ylg=2400;
set(gca,'Visible','off',...
        'ylim',[0 ylg],...
        'xlim',[1993 2016]);


axes('Position',[0.91 0.5 0.08 0.42]);
yy=0;
dy=500;
hold on;
plot([0 0],[0 3000],'-','Color',[0.7 0.3 0]);
for ii=1:6
  y0=(ii-1)*dy+yy;
  plot([-0.2 0.1],[y0 y0],'r-');
  text(0.22,y0,sprintf('%i',y0),'Fontsize',12,'Color','r');
end
set(gca,'xlim',[-0.5 2],...
        'ylim',[0 ylg],...
        'visible','off');

txt{1}=sprintf('Regr: V_{GFWA}= a_0+k*V_{out}');
txt{2}=sprintf('a_0=%5.3g, k=%5.3g',B);
txt{3}=sprintf('95CI a_0: [%5.3g, %5.3g]',bint(1,:));
txt{4}=sprintf('95CI k: [%5.3g, %5.3g]',bint(2,:));
axes('Position',[0.08 0.1 0.5 0.3]);
text(0,0,txt,'Fontsize',12);
set(gca,'visible','off',...
        'ylim',[-0.5 0.5]);

bottom_text(btx,'pwd',1)
  






