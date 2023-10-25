% Plot Extracted T, S, normal U at the specified straits
% Save 2D arrays (depth x width) every N days
%  Corrected flux calculation with T,S,U allocation is used
% see anls_fluxes/extr_TSVdaily_straits08.m
%
% Calcualte fluxes similar to L. de Steur - within 8-2W section of Fram
%

addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/Wavelet;
startup;

close all
clear


f_dataget = 1;

snm = 'FramStr';

lon1=-8.01;
lon2=-1.99;

expt=112;
TV=11;
YR1=2005;
YR2=2019;
dday=7;


pthout = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_straits/';
pthmat = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_theresa/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';

btx='fram_flx_section.m';

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
  jcc=0;
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
%
% Time
    TM = [TM;tmy];


%
% Calculate fluxes for sections in Fram Strait
    if ~exist('il1','var')
      long=SCT(ii0).long;
      dd=abs(long-lon1);
      il1=find(dd==min(dd));
      dd=abs(long-lon2);
      il2=find(dd==min(dd));
    end;

    if ~exist('dZ','var');
      ZMi=SCT(ii0).ZZintrp;
      kz = length(ZMi);
      clear ZZi
      ZZi(1)=0;
      for kk=1:kz-1
        ZZi(kk+1,1)=0.5*(ZMi(kk)+ZMi(kk+1));
      end
      dZ=abs(diff(ZZi));
      dZ(kz)=dZ(kz-1);

%
% Grid cell area along the section
      dL=SCT(ii0).segm_dL;
      Acell=dZ*dL';
    end


    U=SCT(ii0).Unrm;
    T=SCT(ii0).T;
    S=SCT(ii0).S;

    [a1,a2,a3]=size(U);

    for id1=1:a1,
      jcc=jcc+1;
      uu=squeeze(U(id1,:,il1:il2));
      ss=squeeze(S(id1,:,il1:il2));
      tt=squeeze(T(id1,:,il1:il2));
%
% mean S
      II=find(~isnan(ss));
      AA=Acell(:,il1:il2);
      Atot=sum(AA(II));
      smn=nansum(nansum(ss.*AA))./Atot;
      SMN(jcc,1)=smn;

%
% Volume flux:
      vflx=nansum(nansum(uu.*AA));
      VFLX(jcc,1)=vflx;
      
%
% S flux, kT - ~s*1e-3/1 kg of sea water
      rhow=sw_dens0(ss,tt);
      sflx=nansum(nansum((ss*1e-3).*uu.*AA.*rhow))*1e-6; % kT/s
      SFLX(jcc,1)=sflx;
%
% FW flux - overall integration
      srf=34.9;
      fwf1=nansum(nansum( (srf-ss)./srf.*uu.*AA));  % 
      FWF1(jcc,1)=fwf1;
%
% FW flux - only down to sref - ignor S>sref
      ssm=ss;
      ssm(ssm>srf)=nan;
      fwf2=nansum(nansum( (srf-ssm)./srf.*uu.*AA));  % 
      FWF2(jcc,1)=fwf2;


    end;  % days 

  end  % years

end



% Convert to Sv, mSv, kT
% Flux positive south! 
Vflx=-VFLX*1e-6;  % m3/s --> Sv, positive southward
FWflx1=-FWF1*1e-3; % mSv
FWflx2=-FWF2*1e-3; % mSv
Sflx=-SFLX;  % kT/s


% ---------------
% Plotting
% ---------------
%
% Vol Flux
%
sunits='Sv';
tk1=-20;
tk2=20;
dtk=1;
clr1=[0 0.4 0.8];
clr2=[0 0.6 1];
nfg=1;
flxnm='8-2W, VFlux';
sub_plotflx(Vflx, nfg, sunits, tk1, dtk, tk2, TM, clr1, clr2, flxnm);
bottom_text(btx,'pwd',1);


%
% FW flux, Sref=34.9;
% Integrate over all S values
sunits='mSv';
tk1=0;
tk2=250;
dtk=10;
clr1=[0 0.8 0.2];
clr2=[0.4 1 0.6];
nfg=2;
flxnm='8-2W, FWFlux(S), Sref=34.9';
sub_plotflx(FWflx1, nfg, sunits, tk1, dtk, tk2, TM, clr1, clr2, flxnm);
bottom_text(btx,'pwd',1);


%
% FW flux, Sref=34.9;
sunits='mSv';
tk1=0;
tk2=250;
dtk=10;
clr1=[0 0.8 0.2];
clr2=[0.4 1 0.6];
nfg=3;
flxnm='8-2W, FWFlux(zS0), Sref=34.9';
sub_plotflx(FWflx2, nfg, sunits, tk1, dtk, tk2, TM, clr1, clr2, flxnm);
bottom_text(btx,'pwd',1);

%
% South S flux
sunits='kT/s';
tk1=100;
tk2=800;
dtk=50;
clr1=[0.8 0.2 0];
clr2=[1 0.6 0.4];
nfg=4;
flxnm='8-2W, S Flux';
sub_plotflx(Sflx, nfg, sunits, tk1, dtk, tk2, TM, clr1, clr2, flxnm);
bottom_text(btx,'pwd',1);

% Mean S
sunits='psu';
tk1=33;
tk2=39;
dtk=0.05;
clr1=[0.4 0. 0.8];
clr2=[0.8 0.5 1];
nfg=5;
flxnm='Mean S';
sub_plotflx(SMN, nfg, sunits, tk1, dtk, tk2, TM, clr1, clr2, flxnm);
bottom_text(btx,'pwd',1);













