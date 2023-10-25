% Info about expeirments with CICE5 for 0.08 and 0.04 and CICE4 for 0.08
%
% BL99 - old configuration for CICE4 following Bitz & Lipscomb, 1999
%        this includes thermal T profile, conductivity, sh/wave radiative transfer
%        albedo, etc.
%  however, in experiment with CICE5 bubbly conductivity and mushy thermodynamics is used
%  even for BL99 option 
%
% Tmlt/Tfrz - choice of melting and freezing points for calculating ice freeze potential
%    in the ocean-ice coupler
%
% Rheologies: EVP and EAP
%
function EXPT = sub_cice_experiments;

% Old experiment with CICE4 - published with FW tracers and Greenland runoff
ii=1;
EXPT(ii).Nmb      = 110;  % 110 - no Greenland runoff & river clim; 112 - Greenl + monthly runoff
EXPT(ii).cice     = 'CICE4';
EXPT(ii).res      = 0.08;
EXPT(ii).cice_opt = 'BL99Tfrz';
EXPT(ii).dir_outp = '/nexsan/archive/ARCc0.08_110/data/'; % + year_cice/
EXPT(ii).pthmat   = '/Net/kronos/ddmitry/hycom/ARCc0.04/datamat/cice_mnth/';
EXPT(ii).flnm     = '112cice4_outp_tser_BL99Tfrz.mat';
EXPT(ii).rheol    = 'EVP';
EXPT(ii).atm_forc = 'CFSR/CFSv2';
EXPT(ii).day_mean = 0;   % = 1 - mean, =0 - instant.

res = 0.08;
expt = 122;
cice = 'CICE5';
for ik=1:4
  ii=ii+1;
  rheol = 'EVP';
  frc = 'CFSv2';
  if ik==1
    texpt = 'dltEddTfrz';  % delta-Eddington Tmlt=Tfrz in ice melt potential
  elseif ik==2
    texpt = 'BL99Tmlt';   % CCSM3 sh/wave radiation, Tmlt > Tfrz ice heat potent.
  elseif ik==3
    texpt = 'BL99Tfrz';   % CSM3 sh/wave radiation, Tmlt=Tfrz (similar to 0.04 original run)
  elseif ik==4
    texpt = 'dltEddTmlt';
  end

  dir_outp = sprintf('/nexsan/people/ddmitry/hycom/ARCc0.04_022/ARCc0.08_%3.3i_%s',...
              expt,texpt);
  EXPT(ii).Nmb      = expt;
  EXPT(ii).cice     = cice;
  EXPT(ii).res      = res;
  EXPT(ii).cice_opt = texpt;
  EXPT(ii).dir_outp = dir_outp;
  EXPT(ii).pthmat   = '/Net/kronos/ddmitry/hycom/ARCc0.08/datamat/cice_mnth/';
  EXPT(ii).flnm     = sprintf('%3.3icice_outp_tser_%s.mat',expt,texpt);
  EXPT(ii).rheol    = rheol;
  EXPT(ii).atm_forc = frc;
  EXPT(ii).day_mean = 1;   % = 1 - mean, =0 - instant.
end


res = 0.04;
cice = 'CICE5';
for ik=1:4
  ii=ii+1;
  rheol = 'EVP';
  expt = 022;
  frc = 'CFSv2'; % forcing
  if ik==1
    texpt = 'BL99Tfrz';   % Original simulation 
  elseif ik==2
    texpt = 'dltEddTmlt';  % delta-Eddington Tmlt>Tfrz in ice freeze potentia
  elseif ik == 3
    texpt = 'dltEddTmltEAP';
    rheol = 'EAP'; 
  elseif ik == 4
    expt = 023;
    texpt = 'dltEddTmltEAPJRA';
    rheol = 'EAP';
    frc = 'JRA';
  end

  dir_outp = sprintf('/nexsan/people/ddmitry/hycom/ARCc0.04_%3.3i/data/',expt);
  EXPT(ii).Nmb      = expt;
  EXPT(ii).cice     = cice;
  EXPT(ii).res      = res;
  EXPT(ii).cice_opt = texpt;
  EXPT(ii).dir_outp = dir_outp;
  EXPT(ii).pthmat   = '/Net/kronos/ddmitry/hycom/ARCc0.04/datamat/cice_mnth/';
  EXPT(ii).flnm     = sprintf('%3.3icice_outp_tser_%s.mat',expt,texpt);
  EXPT(ii).rheol    = rheol;
  EXPT(ii).atm_forc = frc;
  EXPT(ii).day_mean = 1;   % = 1 - mean, =0 - instant.
end

res   = 0.08;
cice  = 'CICE5';
rheol = 'EAP';
expt  = 123;
frc   = 'JRA';
ii = ii+1;
texpt = 'dltEddTmltEAPJRA';
dir_outp = sprintf('/nexsan/people/ddmitry/hycom/ARCc0.08_%3.3i/data/',expt);
EXPT(ii).Nmb      = expt;
EXPT(ii).cice     = cice;
EXPT(ii).res      = res;
EXPT(ii).cice_opt = texpt;
EXPT(ii).dir_outp = dir_outp;
EXPT(ii).pthmat   = '/Net/kronos/ddmitry/hycom/ARCc0.08/datamat/cice_mnth/';
EXPT(ii).flnm     = sprintf('%3.3icice_outp_tser_%s.mat',expt,texpt);
EXPT(ii).rheol    = rheol;
EXPT(ii).atm_forc = frc;
EXPT(ii).day_mean = 1;   % = 1 - mean, =0 - instant.



% Print all experiments:
fprintf('============  CICE experiments   ===========\n');
nexpt = length(EXPT);
for ii=1:nexpt
  nmexp = EXPT(ii).cice_opt;
  res   = EXPT(ii).res;
  vcice = EXPT(ii).cice;
  rheol = EXPT(ii).rheol;
  frc   = EXPT(ii).atm_forc;
  nmbexp= EXPT(ii).Nmb;

  fprintf('#%2.2i %3.2fHYCOM-%s expt-%3.3i %s %s %s\n',ii,res,vcice,nmbexp,nmexp,rheol,frc);

end
fprintf(' ===========================================  \n');




return
