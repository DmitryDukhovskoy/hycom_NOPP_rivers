% Plot Extracted T, S, normal U at the specified straits
% Save 2D arrays (depth x width) every N days
%  Corrected flux calculation with T,S,U allocation is used
% see anls_fluxes/extr_TSVdaily_straits08.m
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/Wavelet;
startup;

close all
clear


f_dataget = 1;

snm = 'FramStr';

expt=112;
TV=11;
YR1=2005;
YR2=2019;
dday=7;


pthout = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_straits/';
pthmat = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_theresa/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';

btx='plot_VolFWFlux_climat.m';

%tmm=[datenum(2008,1,1):7:datenum(2008,12,31)];
%tmm=tmm(:);
%for isc=1:5
%  SCT(isc).Time=tmm;
%end;


% Process data:
if f_dataget==1
  TM   = [];
  Vflx = [];
  FWflx1 = [];
  FWflx2 = [];
  Hflx1  = [];
  Hflx2  = [];
  for YR = YR1:YR2
    yr=YR;
    dE=datenum(yr,12,31);
    dJ1=datenum(yr,1,1);
    ndays=dE-dJ1+1;
  
%    tmy = [dJ1:dday:dE]';

    fmatout=sprintf('%shycom008_%3.3i_StraitFluxesDay_%4.4i.mat',...
                    pthmat,expt,YR);
    fprintf('Loading %s\n',fmatout);
    load(fmatout);

    nsc = length(SCT);
    if ~exist('ii0','var');
			for ik=1:nsc
				nm=SCT(ik).Name;
				if strncmp(nm,snm,4); ii0=ik; break; end;
			end
    end
%
%  tmy
    tmy = SCT(ii0).Time;
%
% Net Vol flux:
    dmm = SCT(ii0).VolFlx_m3s;
    if length(dmm) ~=length(tmy);
      fprintf('%i Missing dates ??? %i should be %i\n',YR, length(dmm),length(tmy));
      keyboard;
    end
    Vflx = [Vflx;dmm];

%

%Sref1 = 34.8;
% FWflx1: 
    dmm = SCT(ii0).FWflx1_m3s;
    FWflx1 = [FWflx1;dmm];
%
%Sref2 = 34.9;
% FWflx2: 
    dmm = SCT(ii0).FWflx2_m3s;
    FWflx2 = [FWflx2;dmm];
%
%Tref1 = -1.8; % Ref T to calc. H flux
% Hflx1
    dmm = SCT(ii0).Hflx1_W;
    Hflx1 = [Hflx1;dmm];
%
%Tref2 = 0;    % Ref T to calc. H flux
% Hflx2
    dmm = SCT(ii0).Hflx2_W;
    Hflx2 = [Hflx2;dmm];
%
% Time
    TM = [TM;tmy];

  end

end

% Convert to Sv, mSv, TW
Vflx=Vflx*1e-6;  % m3/s --> Sv
FWflx1=FWflx1*1e-3; % mSv
FWflx2=FWflx2*1e-3; % mSv
Hflx1=Hflx1*1e-12;  % TW
Hflx2=Hflx2*1e-12;  % TW


% ---------------
% Plotting
% ---------------
%
% Vol Flux
%
sunits='Sv';
tk1=-10;
tk2=10;
dtk=1;
clr1=[0 0.4 0.8];
clr2=[0 0.6 1];
nfg=1;
flxnm='VFlux';
sub_plotflx(Vflx, nfg, sunits, tk1, dtk, tk2, TM, clr1, clr2, flxnm);
bottom_text(btx,'pwd',1);


%
% FW flux, Sref=34.8;
sunits='mSv';
tk1=-250;
tk2=0;
dtk=10;
clr1=[0 0.8 0.2];
clr2=[0.4 1 0.6];
nfg=2;
flxnm='FWFlux, Sref=34.8';
sub_plotflx(FWflx1, nfg, sunits, tk1, dtk, tk2, TM, clr1, clr2, flxnm);
bottom_text(btx,'pwd',1);


%
% FW flux, Sref=34.9;
sunits='mSv';
tk1=-250;
tk2=0;
dtk=10;
clr1=[0 0.8 0.2];
clr2=[0.4 1 0.6];
nfg=3;
flxnm='FWFlux, Sref=34.9';
sub_plotflx(FWflx2, nfg, sunits, tk1, dtk, tk2, TM, clr1, clr2, flxnm);
bottom_text(btx,'pwd',1);

%
% Heat flux, TW
sunits='TW';
tk1=-240;
tk2=240;
dtk=20;
clr1=[0.8 0.2 0];
clr2=[1 0.6 0.4];
nfg=4;
flxnm='HeatFlx, Tref=-1.8C';
sub_plotflx(Hflx1, nfg, sunits, tk1, dtk, tk2, TM, clr1, clr2, flxnm);
bottom_text(btx,'pwd',1);

% Heat flux, TW
sunits='TW';
tk1=-240;
tk2=240;
dtk=20;
clr1=[0.8 0.2 0];
clr2=[1 0.6 0.4];
nfg=5;
flxnm='HeatFlx, Tref=0C';
sub_plotflx(Hflx2, nfg, sunits, tk1, dtk, tk2, TM, clr1, clr2, flxnm);
bottom_text(btx,'pwd',1);













