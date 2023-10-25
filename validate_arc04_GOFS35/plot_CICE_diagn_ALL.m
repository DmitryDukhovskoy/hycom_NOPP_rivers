% Plot monthly statistics of CICE output
% sea ice melt rate and other CICE output characteristics
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
pEXPT = [2:9];  % experiments to plot - caution with CICE4 - instant. output, not all fields avail

%
% Note: be careful comparing expts 110/112 CICE4 output - instanteneous
% saved at 00hr UZT with daily averaged fields from CICE5 experiments
% Different solar, ln/wave fluxes, ice melt (top) etc. !!!


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


FPLT = [2,4,5,6,7,9,10];  % specify fields to plot

btx = 'plot_CICE_diagn_ALL.m';

%
% Combine all experiments:
icc = 0;
AA = struct;
nexpts = length(EXPT);
for ixx = 1:nexpts  % skip CICE4 expts - only a few fields could be compared
  dmm = pEXPT-ixx;
  if all(dmm); continue; end;

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
  clear TM
  for irc=1:nrc
    fprintf('Record irc=%i\n',irc);
    TM(irc,1) = ICE(irc).dnmb;

		nplt = length(FPLT);
		for ip=1:nplt
      ifld = FPLT(ip);
      fld = FLDS{ifld};
      if isfield(ICE,fld)
        A = sub_iceA(fld,irc,ICE,A);  
      else
        fprintf('expt %i %s Field is missing: %s\n',ixx,texpt,fld);
        stt = sprintf('A.%s_mn(irc,1) = nan;',fld);
        eval(stt)
      end
    end
  end

  icc=icc+1;
  AA(icc).expt_indx = ixx;
  AA(icc).TM = TM;
  AA(icc).A = A;
end

CLR = [1 0.6 0;...
       0.8 0.55 0; ...
       0.5 1   0.2; ...
       0   0.6 0.9; ...
       0.8 0   0.9; ...
       1  0.1  0.1; ...
       0  0.7  0.2; ...
       0.6 0.  0.5; ...
       0  0.   0.9];


dv=datevec(TM(1));
for ip=1:nplt
  ifld = FPLT(ip);
  fld = FLDS{ifld};

  dv1 = datevec(TM(1));
	dv2 = datevec(TM(end));
  sttl = sprintf('%s, %i/%i - %i/%i',fld,dv1(1:2),dv2(1:2));
  fgn = ip;

% Plot all experiments combined:
  figure(ip); clf;
  set(gcf,'Position',[1597 504 919 727]);
  axes('Position',[0.09 0.5 0.83 0.42]);

  Xmin = 1e20;
  Xmax = -1e20;
  Ymin = 1e20;
  Ymax = -1e20;
  
  tm0 = datenum(2017,1,1);
  for icc = 1:length(AA)
    A = AA(icc).A;  
    TM = AA(icc).TM;
    ixx = AA(icc).expt_indx;
    clr = CLR(ixx,:);
    if isempty(A); continue; end

    [ymin,ymax,xmin,xmax] = sub_plot_iceA_all(fld,A,TM,tm0,clr);

    Xmin = min([Xmin,xmin]);
    Xmax = max([Xmax,xmax]);
    Ymin = min([Ymin,ymin]);
    Ymax = max([Ymax,ymax]);

  end

  stl = sprintf('CICE, median, fld=%s',fld);
  title(stl);

  lpos = [0.1 0.1 0.6 0.3];
  TM0 = [tm0+Xmin:tm0+Xmax];
  sub_setgca_iceA_all(ip,EXPT,CLR,TM0,tm0,Ymin,Ymax,lpos);

  bottom_text(btx,'pwd,',1);
end






