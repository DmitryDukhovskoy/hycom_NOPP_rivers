% Plot sea ice melt rate
% extracted in get_daily_melt08.m
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

expt = 022;
pthout = '/Net/kronos/ddmitry/hycom/ARCc0.04/datamat/cice_mnth/';
pthmat = pthout;

%texpt = 'noTmlt'; % in ocean-cice coupler, Tfrz is changed to Tmlt when T>Tmlt
               % to calculate frzmlt sent to CICE 
texpt = 'dltEddTmlt';

FPLT = [1:16];  % specify figures to produce

btx = 'plot_daily_melt04.m';

fmatout = sprintf('%s%3.3icice_outp_tser_%s.mat',pthmat,expt,texpt);
fprintf('Loading %s\n',fmatout);
load(fmatout);

fprintf('Deriving staistics CICE \n');
nrc = length(ICE);
A = struct;
for irc=1:nrc

  dmm = ICE(irc).dnmb;
%  if isempty(dmm);
%    TM(irc,1) = NaN;
%    continue;
%  end

  fprintf('Record irc=%i\n',irc);
  TM(irc,1) = ICE(irc).dnmb;

  A = sub_iceA('ai',irc,ICE,A);      % ice area, fraction
  A = sub_iceA('fswdn',irc,ICE,A);   % down solar flux, W/m2
  A = sub_iceA('frzmlt',irc,ICE,A);  % freeze-melt potential, >0 - ice forms
%  A = sub_iceA('fswint',irc,ICE,A);  % sh/wave absorbed in ice interior, W/m2 <--- all 0s???
  A = sub_iceA('fswabs',irc,ICE,A);  % snow/ice/ocean absorbed solar flux
  A = sub_iceA('fhocn',irc,ICE,A);   % heat flux ice to ocean, W/m2
  A = sub_iceA('fswthru',irc,ICE,A); % SW flux thru ice to ocean
  A = sub_iceA('meltt',irc,ICE,A);   % top ice melt, cm/day
  A = sub_iceA('melts',irc,ICE,A);   % top snow melt, cm/day
  A = sub_iceA('meltb',irc,ICE,A);   % basal ice melt, cm/day
  A = sub_iceA('flwdn',irc,ICE,A);   % down longwave flux, W/m2
  A = sub_iceA('congel',irc,ICE,A);  % congelation ice grwoth, cm/day
  A = sub_iceA('albsni',irc,ICE,A);  % snow/ice broad band albedo
  A = sub_iceA('fsens',irc,ICE,A);   % sensible heat flux, "de-normalized" (not _ai)
  A = sub_iceA('flat',irc,ICE,A);    % latent heat flux, "de-normalized"
  A = sub_iceA('tsfc',irc,ICE,A);    % snow/ice surf T
  A = sub_iceA('sst',irc,ICE,A);     % ocean surface SST
end

dv=datevec(TM(1));

fgn = 1;
if ~all(FPLT-fgn)
		fld = 'fswdn';
		sttl = sprintf('0.04-%3.3i %s %s, %i',expt,texpt,fld,dv(1));
		sub_plot_iceA(fld,A,TM,sttl,fgn);
		bottom_text(btx,'pwd,',1,'Position',[0.04 0.35 0.5 0.05]);
end

fgn = 2;
if ~all(FPLT-fgn)
		fld = 'frzmlt';
		sttl = sprintf('0.04-%3.3i %s %s, %i',expt,texpt,fld,dv(1));
		sub_plot_iceA(fld,A,TM,sttl,fgn);
		bottom_text(btx,'pwd,',1,'Position',[0.04 0.35 0.5 0.05]);
end

fgn = 3;
if ~all(FPLT-fgn)
		fld = 'fswabs';
		sttl = sprintf('0.04-%3.3i %s %s, %i',expt,texpt,fld,dv(1));
		sub_plot_iceA(fld,A,TM,sttl,fgn);
		bottom_text(btx,'pwd,',1,'Position',[0.04 0.35 0.5 0.05]);
end

fgn = 4;
if ~all(FPLT-fgn)
		fld = 'fhocn';
		sttl = sprintf('0.04-%3.3i %s %s, %i',expt,texpt,fld,dv(1));
		sub_plot_iceA(fld,A,TM,sttl,fgn);
		bottom_text(btx,'pwd,',1,'Position',[0.04 0.35 0.5 0.05]);
end

fgn = 5;
if ~all(FPLT-fgn)
		fld = 'fswthru';
		sttl = sprintf('0.04-%3.3i %s %s, %i',expt,texpt,fld,dv(1));
		sub_plot_iceA(fld,A,TM,sttl,fgn);
		bottom_text(btx,'pwd,',1,'Position',[0.04 0.35 0.5 0.05]);
end

fgn = 6;
if ~all(FPLT-fgn)
		fld = 'meltt';
		sttl = sprintf('0.04-%3.3i %s %s, %i',expt,texpt,fld,dv(1));
		sub_plot_iceA(fld,A,TM,sttl,fgn);
		bottom_text(btx,'pwd,',1,'Position',[0.04 0.35 0.5 0.05]);
end

fgn = 7;
if ~all(FPLT-fgn)
		fld = 'melts';
		sttl = sprintf('0.04-%3.3i %s %s, %i',expt,texpt,fld,dv(1));
		sub_plot_iceA(fld,A,TM,sttl,fgn);
		bottom_text(btx,'pwd,',1,'Position',[0.04 0.35 0.5 0.05]);
end

fgn = 8;
if ~all(FPLT-fgn)
		fld = 'meltb';
		sttl = sprintf('0.04-%3.3i %s %s, %i',expt,texpt,fld,dv(1));
		sub_plot_iceA(fld,A,TM,sttl,fgn);
		bottom_text(btx,'pwd,',1,'Position',[0.04 0.35 0.5 0.05]);
end

fgn = 9; 
if ~all(FPLT-fgn)
		fld = 'flwdn';
		sttl = sprintf('0.04-%3.3i %s %s, %i',expt,texpt,fld,dv(1));
		sub_plot_iceA(fld,A,TM,sttl,fgn);
		bottom_text(btx,'pwd,',1,'Position',[0.04 0.35 0.5 0.05]);
end

fgn = 10;
if ~all(FPLT-fgn)
		fld = 'congel';
		sttl = sprintf('0.04-%3.3i %s %s, %i',expt,texpt,fld,dv(1));
		sub_plot_iceA(fld,A,TM,sttl,fgn);
		bottom_text(btx,'pwd,',1,'Position',[0.04 0.35 0.5 0.05]);
end

fgn = 11;
if ~all(FPLT-fgn)
		fld = 'albsni';
		sttl = sprintf('0.04-%3.3i %s %s, %i',expt,texpt,fld,dv(1));
		sub_plot_iceA(fld,A,TM,sttl,fgn);
		bottom_text(btx,'pwd,',1,'Position',[0.04 0.35 0.5 0.05]);
end

fgn = 12;
if ~all(FPLT-fgn)
		fld = 'fsens';
		sttl = sprintf('0.04-%3.3i %s %s, %i',expt,texpt,fld,dv(1));
		sub_plot_iceA(fld,A,TM,sttl,fgn);
		bottom_text(btx,'pwd,',1,'Position',[0.04 0.35 0.5 0.05]);
end

fgn = 13;
if ~all(FPLT-fgn)
		fld = 'flat';
		sttl = sprintf('0.04-%3.3i %s %s, %i',expt,texpt,fld,dv(1));
		sub_plot_iceA(fld,A,TM,sttl,fgn);
		bottom_text(btx,'pwd,',1,'Position',[0.04 0.35 0.5 0.05]);
end

fgn = 14;
if ~all(FPLT-fgn)
		fld = 'tsfc';
		sttl = sprintf('0.04-%3.3i %s %s, %i',expt,texpt,fld,dv(1));
		sub_plot_iceA(fld,A,TM,sttl,fgn);
		bottom_text(btx,'pwd,',1,'Position',[0.04 0.35 0.5 0.05]);
end

fgn = 15;
if ~all(FPLT-fgn)
  fld = 'ai';
  sttl = sprintf('0.04-%3.3i %s %s, %i',expt,texpt,fld,dv(1));
  sub_plot_iceA(fld,A,TM,sttl,fgn);
  bottom_text(btx,'pwd,',1,'Position',[0.04 0.35 0.5 0.05]);
end

fgn = 16;
if ~all(FPLT-fgn)
  fld = 'sst';
  sttl = sprintf('0.04-%3.3i %s %s, %i',expt,texpt,fld,dv(1));
  sub_plot_iceA(fld,A,TM,sttl,fgn);
  bottom_text(btx,'pwd,',1,'Position',[0.04 0.35 0.5 0.05]);
end



