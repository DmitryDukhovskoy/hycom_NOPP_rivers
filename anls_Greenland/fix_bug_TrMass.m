% Extracted monthly Tr Mass in December
% was divided twice by nrec =5 at the end of several year
% due to the bug in extr_MassTr_month.m
%
% Time series of tracer mass by layers
% The Tracer Data are extracted in extr_MassTrcr_mnth.m
% mat files: MassTr01_lrs_*.mat for Greenland tracer

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers
startup;

close all
clear

%FxYR = [1999; 2001; 2002; 2003; 2005; 2006; 2007];
all fixed
nfix = length(FxYR);

nTr = 1; 
nbx = 5; % plot boxes =1,..., nbx

s_fig  = 0;

% Specify levels:
%ilv = 1; % 0-50m
%ilv = 2; % 50-150m
%ilv = 3; % 150-300m
%ilv = 4; % 300-500 m
%ilv = 5; % whole depth 

% Vertical layers
LRS = load('LRS.dat');
nlrs= length(LRS); 


%zz1 = LRS(ilv,1);
%zz2 = LRS(ilv,2);
%
%dz=abs(zz2-zz1);
%
%fprintf('Time Ser. Tracer Mass budget, ilv=%i, %i - %i, nTr=%i\n',ilv,zz1,zz2,nTr);


regn = 'ARCc0.08';
expt = 110;  
pthfig  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/fig_NAtl/',expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
pthsav  = sprintf('/Net/tholia/ddmitry/hycom/ARCc0.08/data_mat%3.3i/',expt);
%fmout   = sprintf('%sGrTrFrct_GrWVol_Regions_NAtl_lev%2.2i.mat',pthmat,ilv);
%btx     = 'propagation_timing.m';

LRS = load('LRS.dat');
nlrs= length(LRS); 

%  clear sumT
for nny = 1:nfix
  iyr = FxYR(nny);
  imo = 12;
  dnmb = datenum(iyr,imo,15);
  dv0  = datevec(dnmb);
  iday = dnmb-datenum(iyr,1,1)+1;

  fmat = sprintf('%sMassTr%2.2i_lrs_%i%2.2i.mat',pthmat,nTr,iyr,imo);
  fprintf('Loading %s\n',fmat);
  if exist(fmat,'file')
    load(fmat);
  else
    fprintf(' =========  MISSING %s\n',fmat);
    error(' Couldnt find file');
  end
    
  for ilv = 1:5
    fprintf('Reading: %i/%2.2i, Tracer %i, Lev %i\n',...
	    iyr,imo,nTr,ilv);
    nrec = TRCR(ilv).nrec;
    Tr = squeeze(TRCR(ilv).MassTr_kg)*nrec;
    TRCR(ilv).MassTr_kg = Tr;
    
    Tr(Tr<1e-12)=nan;
    a=nansum(nansum(Tr));
    fprintf('Lev %i, Corrected tot mass=%10.8d\n',ilv,a)
%    keyboard
  end
  
%  keyboard
  fprintf('Saving corrected %s\n',fmat);
  save(fmat,'TRCR');
end


