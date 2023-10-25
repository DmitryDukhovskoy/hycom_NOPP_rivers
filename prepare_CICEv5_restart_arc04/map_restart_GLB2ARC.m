Code not finished - used Dave Hebert code on gordon
to convert GLBc restart cice5 -> ARCc0.04 cice5

% Map GLBc0.04 GOFS3.5 CICEv5 restart file
% from GLBc --> ARCc0.04
% Possibly need to adjust topo to T17DD
% as GOFS3.5 uses T27

% Read CICE v5 restart

addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
startup

clear all
close

R = 'ARCc0.04';

pthdata = '/Net/kronos/ddmitry/hycom/GOFS3.5/';

fgmapa = sprintf('%sregional.gmapi_GLBc0.04.a',pthdata);
IDM = 3200;
JDM = 5040;
nn  = IDM;
mm  = JDM;

%
% Read gmapi indices
[xmap,ymap] = read_arc_gmapi(IDM,JDM,fgmapa);


fin = sprintf('%scice.restart.2017010109.nc',pthdata);
NC = ncinfo(fin);

%
% Global HYCOM grid/depth
flbase = sprintf('%sregional.grid',pthdata);
fldpth = sprintf('%sdepth_GLBc0.04_27.a',pthdata);
GRD = read_grid_bath(flbase,fldpth);
LON = GRD.PLON;
LAT = GRD.PLAT;


% In GLBc CICE grid is different
% than HYCOM
f_chck_grid=0;
if f_chck_grid==1
% CICE grid: from instanteneous output
  fcice0 = sprintf('%s216_cice5_inst.2017-07-01-00000.nc',pthdata);
  TLON = nc_varget(fcice0,'TLON');
  TLAT = nc_varget(fcice0,'TLAT');

  ih = 1;
  jh = 1;
  ic = 1;
  jc = 1;

  D1 = distance_spheric_coord(TLAT(:,ic),TLON(:,ic),LAT(1:end-1,ih),LON(1:end-1,ih));
  D2 = distance_spheric_coord(TLAT(:,ic),TLON(:,ic),LAT(2:end,ih),LON(2:end,ih));



nvar = length(NC.Variables); 
for iv=1:nvar
  vnm = NC.Variables(iv).Name;
  Agl  = nc_varget(fin,vnm);

  if ~exist('IMAP','var')
% Create mapping index
    [mmg,nng]=size(Ain);
    for ii=1:nn
      JJ=[1:mm]';
      II=ones(mm,1)*ii;
      Iarc = sub2ind([mm,nn],JJ,II);
      Iglb = sub2ind([mmg,nng],ymap(:,ii),xmap(:,ii));

        ja = sub2ind([mm,nn],

end





