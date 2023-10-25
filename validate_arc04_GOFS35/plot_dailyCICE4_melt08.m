% CICE4 output fields
% for GOFS3.1 HYCOM-CICE experiments (11.2) 
% 
% Plot sea ice melt rate
% extracted in get_dailyCICE4_melt08.m
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

expt = 112;
pthout = '/Net/kronos/ddmitry/hycom/ARCc0.04/datamat/cice_mnth/';
pthmat = pthout;

texpt = 'noTmlt'; % Tmlt in ocean-cice coupler, Tfrz is changed to Tmlt when T>Tmlt
               % to calculate frzmlt sent to CICE5

btx = 'plot_dailyCICE4_melt08.m';

fmatout = sprintf('%s%3.3icice_outp_tser_%s.mat',pthmat,expt,texpt);
fprintf('Loading %s\n',fmatout);
load(fmatout);

fprintf('Deriving staistics CICE \n');
nrc = length(ICE);
A = struct;
for irc=1:nrc
  TM(irc,1) = ICE(irc).dnmb;

  if isempty(ICE(irc).ai)
    A.ai_mn(irc,1) = nan;
  end

  A = sub_iceA('ai',irc,ICE,A);      % ice area, fraction
  A = sub_iceA('fswdn',irc,ICE,A);   % down solar flux, W/m2
  A = sub_iceA('frzmlt',irc,ICE,A);  % freeze-melt potential, >0 - ice forms
%  A = sub_iceA('fswint',irc,ICE,A);  % sh/wave absorbed in ice interior, W/m2 <--- all 0s???
%  A = sub_iceA('fswabs',irc,ICE,A);  % snow/ice/ocean absorbed solar flux
%  A = sub_iceA('fhocn',irc,ICE,A);   % heat flux ice to ocean, W/m2
%  A = sub_iceA('fswthru',irc,ICE,A); % SW flux thru ice to ocean
  A = sub_iceA('meltt',irc,ICE,A);   % top ice melt, cm/day
%  A = sub_iceA('melts',irc,ICE,A);   % top snow melt, cm/day
  A = sub_iceA('meltb',irc,ICE,A);   % basal ice melt, cm/day
  A = sub_iceA('flwdn',irc,ICE,A);   % down longwave flux, W/m2
  A = sub_iceA('congel',irc,ICE,A);  % congelation ice grwoth, cm/day
%  A = sub_iceA('albsni',irc,ICE,A);  % snow/ice broad band albedo
  A = sub_iceA('fsens',irc,ICE,A);   % sensible heat flux, "de-normalized" (not _ai)
  A = sub_iceA('flat',irc,ICE,A);    % latent heat flux, "de-normalized"
  A = sub_iceA('tsfc',irc,ICE,A);    % snow/ice surf T

end

dv=datevec(TM(1));

fgn = 1;
fld = 'fswdn';
sttl = sprintf('0.08-%3.3i %s %s, %i',expt,texpt,fld,dv(1));
sub_plot_iceA(fld,A,TM,sttl,fgn);
bottom_text(btx,'pwd,',1,'Position',[0.04 0.35 0.5 0.05]);

fgn = 2;
fld = 'frzmlt';
sttl = sprintf('0.08-%3.3i %s %s, %i',expt,texpt,fld,dv(1));
sub_plot_iceA(fld,A,TM,sttl,fgn);
bottom_text(btx,'pwd,',1,'Position',[0.04 0.35 0.5 0.05]);

fgn = 5;
fld = 'meltb';
sttl = sprintf('0.08-%3.3i %s %s, %i',expt,texpt,fld,dv(1));
sub_plot_iceA(fld,A,TM,sttl,fgn);
bottom_text(btx,'pwd,',1,'Position',[0.04 0.35 0.5 0.05]);

fgn = 6;
fld = 'meltt';
sttl = sprintf('0.08-%3.3i %s %s, %i',expt,texpt,fld,dv(1));
sub_plot_iceA(fld,A,TM,sttl,fgn);
bottom_text(btx,'pwd,',1,'Position',[0.04 0.35 0.5 0.05]);

fgn = 8; 
fld = 'flwdn';
sttl = sprintf('0.08-%3.3i %s %s, %i',expt,texpt,fld,dv(1));
sub_plot_iceA(fld,A,TM,sttl,fgn);
bottom_text(btx,'pwd,',1,'Position',[0.04 0.35 0.5 0.05]);

fgn = 9;
fld = 'congel';
sttl = sprintf('0.08-%3.3i %s %s, %i',expt,texpt,fld,dv(1));
sub_plot_iceA(fld,A,TM,sttl,fgn);
bottom_text(btx,'pwd,',1,'Position',[0.04 0.35 0.5 0.05]);

fgn = 11;
fld = 'fsens';
sttl = sprintf('0.08-%3.3i %s %s, %i',expt,texpt,fld,dv(1));
sub_plot_iceA(fld,A,TM,sttl,fgn);
bottom_text(btx,'pwd,',1,'Position',[0.04 0.35 0.5 0.05]);

fgn = 12;
fld = 'flat';
sttl = sprintf('0.08-%3.3i %s %s, %i',expt,texpt,fld,dv(1));
sub_plot_iceA(fld,A,TM,sttl,fgn);
bottom_text(btx,'pwd,',1,'Position',[0.04 0.35 0.5 0.05]);

fgn = 13;
fld = 'tsfc';
sttl = sprintf('0.08-%3.3i %s %s, %i',expt,texpt,fld,dv(1));
sub_plot_iceA(fld,A,TM,sttl,fgn);
bottom_text(btx,'pwd,',1,'Position',[0.04 0.35 0.5 0.05]);

fgn = 15;
fld = 'ai';
sttl = sprintf('0.08-%3.3i %s %s, %i',expt,texpt,fld,dv(1));
sub_plot_iceA(fld,A,TM,sttl,fgn);
bottom_text(btx,'pwd,',1,'Position',[0.04 0.35 0.5 0.05]);




