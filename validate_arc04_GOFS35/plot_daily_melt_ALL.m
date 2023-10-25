% Plot sea ice melt rate and other CICE output characteristics
% extracted in get_daily_meltXXX.m
% For different experiemnts and CICE4/CICE5 and 0.08 and 0.04
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

EXPT = sub_cice_experiments; 


FLDS{1} = 'ai'; % ice area, fraction
FLDS{2} = 'fswdn';   % down solar flux, W/m2
FLDS{3} = 'frzmlt';  % freeze-melt potential, >0 - ice forms
FLDS{4} = 'fswabs';  % snow/ice/ocean absorbed solar flux
FLDS{5} = 'fhocn';   % heat flux ice to ocean, W/m2
FLDS{6} = 'fswthru'; % SW flux thru ice to ocean
FLDS{7} = 'meltt';   % top ice melt, cm/day
FLDS{8} = 'melts';   % top snow melt, cm/day
FLDS{9} = 'meltb';   % basal ice melt, cm/day
FLDS{10} = 'flwdn';   % down longwave flux, W/m2
FLDS{11} = 'congel';  % congelation ice grwoth, cm/day
FLDS{12} = 'albsni';  % snow/ice broad band albedo
FLDS{13} = 'fsens';   % sensible heat flux, "de-normalized" (not _ai)
FLDS{14} = 'flat';    % latent heat flux, "de-normalized"
FLDS{15} = 'tsfc';    % snow/ice surf T
FLDS{16} = 'sst';     % ocean surface SST


FPLT = [9];  % specify fields to plot

btx = 'plot_daily_melt_ALL.m';

%
% Combine all experiments:
AA = struct;
nexpts = length(EXPT);
for ixx = 1:nexpts
  pthmat = EXPT(ixx).pthmat;
  expt   = EXPT(ixx).Nmb;
  texpt  = EXPT(ixx).cice_opt;
  flnm   = EXPT(ixx).flnm;

  fmatout = sprintf('%s%s',pthmat,flnm);
  fprintf('Loading %s\n',fmatout);
  load(fmatout);

  fprintf('Deriving staistics CICE \n');
  nrc = length(ICE);
  A = struct;
  for irc=1:nrc
    fprintf('Record irc=%i\n',irc);
    TM(irc,1) = ICE(irc).dnmb;

		nplt = length(FPLT);
		for ip=1:nplt
      ifld = FPLT(ip);
      fld = FLDS{ifld};
      A = sub_iceA(fld,irc,ICE,A);  
    end
  end

  AA(ixx).TM = TM;
  AA(ixx).A = A;
end


dv=datevec(TM(1));
for ip=1:nplt
  ifld = FPLT(ip);
  fld = FLDS{ifld};

  dv1 = datevec(TM(1));
	dv2 = datevec(TM(end));
  sttl = sprintf('%s, %i/%i - %i/%i',fld,dv1(1:2),dv2(1:2));
  fgn = ip;

  for ixx=1:nexpts
		A = AA(ixx).A;  
    TM = AA(ixx).TM;
		fstart = ixx;
    if ixx == nexpts, fstart = -1; end;
    sub_plot_iceA_all(fld,A,TM,sttl,fgn,fstart);
  end

  bottom_text(btx,'pwd,',1,'Position',[0.04 0.35 0.5 0.05]);

end




