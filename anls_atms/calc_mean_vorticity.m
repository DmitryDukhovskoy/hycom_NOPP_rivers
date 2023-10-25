% Calculate area-mean vorticity
% using circulation theorem
% similar to Bourassa
%
% from NCEP CFSR and CFSv2
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

RR   = 50; % ring size (diamter), km  
%regn = 'natl'; % Natl region
%regn = 'arctA'; % 
regn = 'arctB'; % 
%regn = 'arctic'; % 

s_mat = 1;

pthdat1 = '/Net/kronos/ddmitry/ncep_cfsr/';
pthdat2 = '/Net/kronos/ddmitry/ncep_cfsv2/';
pthmat  = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.72/data_mat/';
pth72   = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.72/topo_grid/';
pthout  = ('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/data_mat/');
%fmat    = sprintf('%smonthly_areamean_vort_%s%3.3ikm.mat',pthout,regn,RR);
fmat    = sprintf('%smonthly_areamean_vort_%s%3.3ikm.mat',pthout,regn,RR);

btx = 'calc_mean_vorticity.m';

switch(regn),
 case('natl');
%  ii1=41;
%  ii2=122;
%  jj1=9;
%  jj2=113;
  IJr = [41   9
	 41  113
	 122 113
	 41   9];
  hmin= -200;
 case('arctA');
  IJr=[ 67   212
    75   213
   130   198
   167   133
   167    87
   142    78
   122   103
    96   108
    81   113
    30   169];
  hmin = -10;
 case ('arctB');
  IJr = [    52   185
    69   189
    90   187
   118   174
   132   169
   131   136
   121   110
   106   108
    87   122
    58   156
    46   170
    48   181];
  hmin=-10;
 case('arctic'); % exclude near-coastal regions
  IJr=[     61   188
    75   187
    90   189
   107   175
   124   175
   130   174
   134   167
   132   154
   130   149
   133   144
   135   139
   143   138
   145   128
   144   121
   148   113
   151   107
   152    94
   147    89
   139    85
   122   104
   103   110
    96   118
    78   130
    65   141
    59   151
    51   163
    47   173
    48   180
    56   186]; 
  hmin = -10;
end
  

% Get indices for interpolating onto HYCOM ARCc0.72 grid:
fgrd = sprintf('%scfsr_gridindx_arc072_nghb.mat',pthmat);
load(fgrd);
mc = INDX.dim_rows;
nc = INDX.dim_colms;

% Topo for ARCc0.72
fsv=[pth72,'new_bath072.mat'];
load(fsv);
LN=elon;
LT=alat;
HH=hnew;
HH(isnan(HH))=10;
clear hnew elon alat;
[mm,nn]=size(HH);
[DX,DY]=sub_dx_dy(LN,LT);
[II,JJ] = meshgrid([1:nn],[1:mm]);

% Define points in the region
inp = inpolygon(II,JJ,IJr(:,1),IJr(:,2));

% Define deep ocean points only
Ioc = find(inp==1 & HH<hmin);
Inc = find(inp~=1 & HH>=hmin);


% Get PANG - angles to rotate vectors into HYCOM grid
% b) Rotate the velocities HYCOM ARCc to x-wards y-wards.  
%The array pang in regional.grid can be used to do this.  
%The eastwards,northwards to x-wards,y-wards code is in ALL/force/src/wi.f:
%
%              COSPANG  = COS(PANG(I,J))
%              SINPANG  = SIN(PANG(I,J))
%              TXM(I,J) = COSPANG*TXMIJ + SINPANG*TYMIJ
%              TYM(I,J) = COSPANG*TYMIJ - SINPANG*TXMIJ
fina = sprintf('%sregional.grid.a',pth72);
finb = sprintf('%sregional.grid.b',pth72);
[PANG,nn,mm] = read_pang(fina,finb);
COSPANG = cos(PANG);
SINPANG = sin(PANG);

cc=0;
CRL=[];
CINDX = [];

VRTC.Title='Area Mean Vorticity Ocean Points Monthly';
VRTC.Info = 'CFSR/CFS fields interpolated onto ARCc0.72 grid';
VRTC.Index_Ocean=Ioc;
VRTC.Size_Domain=[mm,nn];
VRTC.File_Grid=fgrd;
VRTC.Ring_Size=RR;

for yr=1993:2016
  
  if yr>=2011
    fgrd = sprintf('%scfs_gridindx_arc072_nghb.mat',pthmat);
    load(fgrd);
    mc = INDX.dim_rows;
    nc = INDX.dim_colms;
  end
  
  if yr<2011
    fnm = sprintf('%scfsr-sea_%i_mon_uv-10m.nc',pthdat1,yr);
  else
    fnm = sprintf('%scfsv2-sec2_%i_mon_uv-10m.nc',pthdat2,yr);
  end

  if ~exist('X','var')
    X  = nc_varget(fnm,'Longitude');
    Y  = nc_varget(fnm,'Latitude');
    Ix = find(X>180);
    X(Ix)=X(Ix)-360;
    [XX,YY]=meshgrid(X,Y);
  end
  dm = nc_varget(fnm,'MT');
  TM=datenum(1900,12,31)+dm;

% Northward and eastward wind components
  U = nc_varget(fnm,'wndewd');
  V = nc_varget(fnm,'wndnwd');
  
  for im=1:12
    tic
    cc = cc+1;
    dnmb = TM(im);
    DtV = datevec(dnmb);
    fprintf('%i/%2.2i/%2.2i\n',DtV(1:3));
    
    uu=squeeze(U(im,:,:));
    vv=squeeze(V(im,:,:));
    s=sqrt(uu.^2+vv.^2);

    uh = uu(INDX.II);
    vh = vv(INDX.II);
    ur = COSPANG.*uh+SINPANG.*vh;
    vr = COSPANG.*vh-SINPANG.*uh;
    sr = sqrt(ur.^2+vr.^2);
    
% Mean curl over Nordic Seas getpts
    [VORT,CINDX] = vort_area_mean(XX,YY,ur,vr,RR,CINDX);

    VRT = VORT.VRT;
    
% Keep only ocean pnts within the region of interest
    dmm = VORT.VRT(Ioc);
    VRTC.VRT(cc,:) = dmm;
    VRTC.TM(cc,1)  = dnmb;

    fprintf('Processing 1 time rec: %8.3f sec\n',toc); 
    
  end
  
  if s_mat==1
    fprintf('======   Saving %s\n',fmat);
    save(fmat,'VRTC');
  end
  
  
  
end


