% Save lon/lat and inidices of the sections
% for python analysis
% obsloc_Fram_FWyr.py
%
% Extracted T, S, normal U at the specified straits
% for processing in python
%  Corrected flux calculation with T,S,U allocation is used
% see anls_fluxes/extr_TSVdaily_straits08.m
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/Wavelet;
startup;

close all
clear

f_save = 1;
snm = 'FramStr';

expt=112;
TV=11;
YR1=2005;
%YR2=2019;
dday=7;


pthout = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_straits/';
pthmat = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_theresa/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';


% Process data:
TM   = [];
Vflx = [];
FWflx1 = [];
FWflx2 = [];
Hflx1  = [];
Hflx2  = [];
YR = YR1;
yr=YR;
dE=datenum(yr,12,31);
dJ1=datenum(yr,1,1);
ndays=dE-dJ1+1;

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
tmy = SCT(ii0).Time;
%
dL = SCT(ii0).segm_dL; 
Dist = cumsum(dL);
lonsct = SCT(ii0).long;
latsct = SCT(ii0).latd;
Indx   = SCT(ii0).IJ_indx;
Jndx   = SCT(ii0).IJ_indx;

nm = SCT(ii0).Name;
fflux = sprintf('%sarc08_sect_coord_%s.dat',pthout,nm);
fprintf('Writing coordinates --> %s\n',fflux);
fid = fopen(fflux,'w','ieee-be');
nfields = 5;
npnts = length(Indx);
fwrite(fid,nfields,'int');
fwrite(fid,npnts,'int');
fwrite(fid,dL,'float32');
fwrite(fid,Indx,'float32');
fwrite(fid,Jndx,'float32');
fwrite(fid,lonsct,'float32');
fwrite(fid,latsct,'float32');
fclose(fid);








