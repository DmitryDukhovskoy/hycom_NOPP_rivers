% Plot extracted average T/S profiles
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

% Specify field to plot:
pfld  = 'temp';
%pfld  = 'salin';

% Specify grid resolution to plot:
ires = 0.04;
%ires = 0.08;


regn = sprintf('ARCc%3.2f',ires);
if ires == 0.08
  expt = 123;
else
  expt = 23;
end
expt_nm = 'dltEddTmltEAPJRA';

hg    = 2^100;
rg    = 9806;
% Change Tv and nlev for different ARCc fields
Tv    = 11; % bathym v11
nlev  = 41; % 41

pthout = '/nexsan/people/ddmitry/hycom/data_extracted/';

fmat_out = sprintf('%shycom%3.3i_%3.3i_meanprof_%s.mat',...
                    pthout,ires*100,expt,pfld);

load(fmat_out);
TM = PROF.TM;
PrfEU = PROF.prof_euras;
PrfCA = PROF.prof_canad;
ZZi = PROF.Depths;

figure('Position',[1639         491         864         850]);
axes('Position',[0.1, 0.1, 0.35, 0.8]);
hold on
plot(PrfEU,ZZi)
set(gca,'xgrid','on',...
        'ygrid','on');
title(sprintf('%4.2f Mean %s Euras',ires,pfld));

axes('Position',[0.55, 0.1, 0.35, 0.8]);
hold on
plot(PrfCA,ZZi)
set(gca,'xgrid','on',...
        'ygrid','on');
title(sprintf('%4.2f Mean %s Canad',ires,pfld));




