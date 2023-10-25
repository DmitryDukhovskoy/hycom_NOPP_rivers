% Calculate upwelling index
% along Greenland coast
% using Bakun method
% M = 1/(f*rho_w)*(tau x k)
% see also Picket and Padun,
% upwelling Cal. Current, JGR
%
% from NCEP CFSR and CFSv2
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

rhoa = 1.2;
Cd   = 0.0013;
rhow = 1025;

s_mat = 1;

pthdat1 = '/Net/kronos/ddmitry/ncep_cfsr/';
pthdat2 = '/Net/kronos/ddmitry/ncep_cfsv2/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/data_mat/';

pthout  = ('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/data_mat/');
fmat    = sprintf('%sarc008_Greenl_upwl_CFSR_month.mat',pthout);

% Get indices for interpolating onto HYCOM ARCc0.72 grid:
fgrd = sprintf('%scfsr_gridindx_arc008_nghb.mat',pthmat);
load(fgrd);
mc = INDX.dim_rows;
nc = INDX.dim_colms;

% Topo for ARCc0.08
ftopo = sprintf('%s/depth_ARCc0.08_11.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
LN = nc_varget(ftopo,'Longitude');
LT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(HH);
[DX,DY]=sub_dx_dy(LN,LT);
[II,JJ] = meshgrid([1:nn],[1:mm]);
Fcor = 2*7.292e-5*sind(LT);

% Load Greenland coast
% with local normal unit vectors
GRC = sub_greenl_coast08(HH);
ngr = length(GRC.X);


% Get PANG - angles to rotate vectors into HYCOM grid
% b) Rotate the velocities HYCOM ARCc to x-wards y-wards.  
%The array pang in regional.grid can be used to do this.  
%The eastwards,northwards to x-wards,y-wards code is in ALL/force/src/wi.f:
%
%              COSPANG  = COS(PANG(I,J))
%              SINPANG  = SIN(PANG(I,J))
%              TXM(I,J) = COSPANG*TXMIJ + SINPANG*TYMIJ
%              TYM(I,J) = COSPANG*TYMIJ - SINPANG*TXMIJ
fina = sprintf('%sregional.grid.a',pthtopo);
finb = sprintf('%sregional.grid.b',pthtopo);
[PANG,nn,mm] = read_pang(fina,finb);
COSPANG = cos(PANG);
SINPANG = sin(PANG);

cc=0;
UPW.Title='Upwelling Index, Bakun, Greenland Coast';
UPW.Info = 'CFSR/CFS fields interpolated onto ARCc0.72 grid';
UPW.Size_Domain=[mm,nn];
UPW.File_Grid=fgrd;

for yr=1993:2016
  if yr>=2011
    fgrd = sprintf('%scfs_gridindx_arc008_nghb.mat',pthmat);
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
%    Fcor = 2*7.292e-5*sind(YY);
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
    
% Upwelling index, Ekman transport
% m3/s
    tx = rhoa*Cd*sr.*ur;
    ty = rhoa*Cd*sr.*vr;
    Mx = 1./(rhow*Fcor).*ty;
    My = -1./(rhow.*Fcor).*tx;
    
% Project on Greenland coast points:
   Mnrm = 0;
   for ik=1:ngr
     i0   = GRC.Iocn(ik);
     j0   = GRC.Jocn(ik);
     xnrm = GRC.xnrm_offshore(ik);
     ynrm = GRC.ynrm_offshore(ik);
     mx   = Mx(j0,i0);
     my   = My(j0,i0);
     mnrm = mx*xnrm+my*ynrm;
     Mnrm(ik) = mnrm;
   end

   UPWL.TM(cc) = dnmb;
   UPWL.Normal_m3s(cc,:) = Mnrm;
   
   f_plt=0;
   if f_plt==1
     fprintf('Plotting ...\n');

     fn = 2;
     sub_plot_upwelling(HH,LN,LT,fn,Mnrm,GRC,dnmb)
     
     btx = 'calc_upwelling.m';
     bottom_text(btx,'pwd',1);      
   end

  end
end

if s_mat==1
  fprintf('Saving %s\n',fmat);
  save(fmat,'GRC','UPWL');
end



