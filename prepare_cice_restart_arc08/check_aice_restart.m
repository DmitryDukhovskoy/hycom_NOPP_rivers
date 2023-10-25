% Reads restart for CICE v5 
% netCDF file
% 
% ice_itd.F90 check aggregated aice - sends an error if aice > 1
%
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close

YH = 117;
AB = 'f';

ncat = 5;

pthin = '/nexsan/people/ddmitry/hycom/ARCc0.04_022/incoming/';
%restart created map_arc04_to_arc08.m:
finp = sprintf('%scice.restart08_%3.3i%s_DD.nc',pthin,YH,AB);
%finp = sprintf('%scice.restart04_%3.3i%s.nc',pthin,YH,AB);

% Ice concentration by cateogries:
aicen = nc_varget(finp,'aicen');
% Total
Ai = squeeze(sum(aicen,1));

% Check aggregate ice area:
eps=1e-22;
I=find(Ai>1.0);

% Correct excess ice:
% Implemented in CICE: ice_itd.F90
f_crct=0;
if f_crct==1
  if ~isempty(I)
  % adjust aice by categories:
    fprintf('Found %i out of bound aicen\n',length(I));
    for ii=1:length(I);
      i0=I(ii);
      dlt = Ai(i0)-1+eps;
      ww=squeeze(aicen(:,i0)/Ai(i0));
      crc=dlt*ww;
      fprintf('Excess %9.6g, aice=%8.6g\n',dlt,Ai(i0));
      for in=1:ncat
        aicen(in,i0)=aicen(in,i0)-crc(in);
      end
      atot = squeeze(sum(aicen(:,i0)));
      dlt = atot-1;
      fprintf('Corrected: Excess %9.6g, aice=%8.6g\n',dlt,atot);
    end
  end
end

