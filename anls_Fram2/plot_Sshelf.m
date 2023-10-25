% Plot FW content shelf vs deep Fram
%
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
  FWC = [];
  cc=0;
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
%
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

    S = SCT(ii0).S;
    [a1,a2,a3]=size(S);
    Srf1 = 34.9;
    for jj=1:a1
      dmm=squeeze(S(jj,:,:));
      dmm(dmm>Srf1)=nan;
      fwc = nansum((Srf1-dmm)/Srf1.*Acell); % m3 of FW
%      fwc=fwc(:);
%      fwc = fwc./dL;    % m of FW
      cc=cc+1;
      FWC(cc,:)=fwc;
    end
%
% Time
    TM = [TM;tmy];

  end

end


% ---------------
% Plotting
% ---------------
%
% FWC along the section
%

xlon=SCT(ii0).long;
dL=SCT(ii0).segm_dL;


DV = datevec(TM);
Ism=find(DV(:,2)>5 & DV(:,2)<=10);
Iwnt=find(DV(:,2)<=5 | DV(:,2)>10);

dmm=FWC(Ism,:);

fws=nanmean(dmm); 
fws=fws(:);
cumfws=cumsum(fws);  % cum FWC, m3
fws=fws./dL;  % m of FW
ps1=prctile(dmm,25);
ps1=ps1(:);
ps1=ps1./dL;
ps2=prctile(dmm,75);
ps2=ps2(:);
ps2=ps2./dL;

% Winter
dmm=FWC(Iwnt,:);

fww=nanmean(dmm); 
fww=fww(:);
cumfww=cumsum(fww);  % cum FWC, m3
fww=fww./dL;  % m of FW
pw1=prctile(dmm,25);
pw1=pw1(:);
pw1=pw1./dL;
pw2=prctile(dmm,75);
pw2=pw2(:);
pw2=pw2./dL;


figure(1); clf;
axes('Position',[0.08 0.6 0.85 0.32]);
hold on;
plot(xlon,fww,'-','Color',[0 0.2 0.6],'Linewidth',3);
plot(xlon,fws,'-','Color',[0 0.8 0.2],'Linewidth',3);
lgd = legend('NDJFMAM','JJAS');
set(lgd,'Position',[0.7 0.8 0.15 0.09],'Fontsize',14);
%plot(xlon,pw1);
%plot(xlon,pw2);

title('0.08 HYCOM-CICE-110, FWC (Sref=34.9), m, Fram Section');
set(gca,'xlim',[-20 12],...
        'tickdir','out',...
        'ylim',[0 15.5],...
        'xtick',[-20:2:20],...
        'ytick',[0:2:20],...
        'xgrid','on',...
        'ygrid','on',...
        'Fontsize',14);


axes('Position',[0.08 0.1 0.85 0.32]);
hold on
plot(xlon,cumfww,'-','Color',[0 0.2 0.6],'Linewidth',3);
plot(xlon,cumfws,'-','Color',[0 0.8 0.2],'Linewidth',3);

title('intgr{ FWC (Sref=34.9)}, m3, Fram Section'); % FWC m3 per 1 m width section = m2 of FW
set(gca,'xlim',[-20 12],...
        'tickdir','out',...
        'ylim',[0 5e6],...
        'xtick',[-20:2:20],...
        'ytick',[0:1e6:5e6],...
        'xgrid','on',...
        'ygrid','on',...
        'Fontsize',14);

btx = 'plot_Sshelf.m';
bottom_text(btx,'pwd',1);













