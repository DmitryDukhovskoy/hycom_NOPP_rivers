% Calculate tracer fluxes and tracer budget
% at the OB
% to analyze Tracer fluxes
% from the OB into the domain
%
% mat files are saved by years
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

yr1 = 2013;
yr2 = 2015;

s_mat = 0;% =2 - load saved and start from the last record, save mat
          % =1 - start from day 1 - end, save mat
	  % =0 - start from day 1 - end, do not dump mat file

if s_mat==0,
  fprintf('Mat file is not created\n');
elseif s_mat == 1
  fprintf('Mat file will be saved, old mat file will be overridden\n');
elseif s_mat == 2
  fprintf('Extraction is skipped for months where old mat files exist\n');
end

fprintf('======    Years: %i-%i\n',yr1,yr2);

rg = 9806; 
%cBr = 14; % coefficient for Bering Str. tracer - not right
cBr = sub_BeringTrConc(0);% coefficient for Bering Str., tracer, kg/m3

regn = 'ARCc0.08';
expt = 110;  
%pthfig  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/fig_trac/',expt);
pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);


ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2

hmsk=HH;
hmsk(HH<0)=nan;

% Specify segments of the x-section
i1=193;
i2=1021;
j1=21;
j2=j1;
IJs=[ i1  j1;
      i2  j2];
Rmsk = HH;
Rmsk(j1+1:end,:)=0;
Rmsk(:,1:i1)=0;
Rmsk(:,i2+1:end)=0;

nij=size(IJs,1);
IIs=[];
JJs=[];
for ii=1:nij-1
  i1=IJs(ii,1);
  i2=IJs(ii+1,1);
  j1=IJs(ii,2);
  j2=IJs(ii+1,2);
  [I,J]=sub_xsct_indx(i1,j1,i2,j2);
  if size(I,1)==1;
    I=I';
    J=J';
  end
  
  IIs=[IIs;I];
  JJs=[JJs;J];
end;

IJs=[IIs,JJs];

nS=length(IIs);
clear Xl Yl
for ii=1:nS
  i0=IJs(ii,1);
  j0=IJs(ii,2);
  Xl(ii,1)=LON(j0,i0);
  Yl(ii,1)=LAT(j0,i0);
end;
INDs=sub2ind(size(HH),JJs,IIs);

SEGM.I = IJs(:,1);
SEGM.J = IJs(:,2);
SEGM.Ind = INDs;
SEGM.dx  = DX(INDs);
SEGM.nx  = 0; % norm, x component
SEGM.ny  = 1; % norm, y component


f_map=0;
if f_map>0
  figure(1); clf;
  hold on
  contour(HH,[0 0],'k');
  contour(HH,[-5000:1000:-100],'Color',[0.7 0.7 0.7]);
  plot(IIs,JJs,'b.-');
end

ip1=1;
mold = 0;
dday = 7; 
for iyr = yr1:yr2
  TRCR = [];
  TRCR = struct;
  yr=iyr;
  fmat = sprintf('%strcrFlx_AtlOB_daily_%4.4i.mat',pthmat,iyr);
  id1=1;
  nrec = 0;
  for iday = id1:dday:365
    tic
    pthbin = sprintf('/nexsan/archive/%s_%3.3i/data/%i/',regn,expt,yr)

    dnmb=datenum(yr,1,1)+iday-1;
    DV=datevec(dnmb);
    imo=DV(2);
    fprintf('Processing Month %i\n',imo);

    fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
    finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);
  
    if ~exist(fina,'file');
      fprintf('Not found: %s\n\n',fina);
      continue;
    end
    

    nrec = nrec+1;
    nTr = 1; % Greenland only

    tic; 
    
    fprintf('Tr# %i, %4.4i_%2.2i_%2.2i: %s\n',...
	nTr,DV(1:3),fina);

    [F,nn,mm,ll] = read_hycom(fina,finb,'tracer','r_tracer',nTr);
    F(F>1e6)=nan;
    Tr = F;

%    [F,n,m,l] = read_hycom(fina,finb,'u-vel.');
%    F(F>1e6)=0;
%    UN=squeeze(F);
%    
    [F,n,m,l] = read_hycom(fina,finb,'v-vel.');
    F(F>1e6)=0;
    VN=squeeze(F);

% Check collocation of variables:
% but land values in U and V should be nan
%    colloc = check_collocatedU(HH,UN,VN);    

    [F,n,m,l] = read_hycom(fina,finb,'thknss');
    F=F./rg;
    F(F>1e10)=0;
    dH=squeeze(F); 
    nlr=l;
    
    
% Fluxes are calculated as:
%
%     |             |
%     |      i,j    |
%     -U(i,j) *     |
%     |       Tr    |
%     |             |
% -------- | V(i,j)-------
% Calculate transport
% Sign convention:
% U+ - positive X (right)
% V+ - positive Y (up)
%keyboard
    I  = SEGM.I;
    J  = SEGM.J;
    DX = SEGM.dx;
    nx = SEGM.nx;
    ny = SEGM.ny;

    ll=length(I);
    clear FLX
    FLX=zeros(nlr,ll-1);
%    keyboard
%    parfor l=1:ll-1   % little segments of a section
    for l=1:ll-1   % little segments of a section
      i1=I(l);
      i2=I(l+1);
      j1=J(l);
      j2=J(l+1);
      dd=DX(l);  % segment length, m
% Find Normal:
% Sign convention:
      di=abs(i2-i1);
      dj=abs(j2-j1);
      if di>0 & dj>0
	error('Not step-like section ...');
      end

% Check norms:
      if (I(1)-I(end))==0 & ny~=0
	error('Check norms: ny not 0, di=0');
      elseif (J(1)-J(end))==0 & nx~=0
	error('Check norms: nx not 0, dj=0');
      end
      
% Check land segments:      
      H1  = HH(j1,i1);
      H2  = HH(j2,i2);
      flx = zeros(nlr,1); 
      U0  = [];
      V0  = [];
      if H1==0 & H2==0; FLX(:,l)=flx; continue; end;
%  Collocate U,V, tracers, 
% interpolate into the middle of the grid
      if di==0  % Y-section, V*norm=0
	Trj1  = squeeze(Tr(:,j1,i1-1));
	Trj2  = squeeze(Tr(:,j1,i1));
	U0    = squeeze(UN(:,j1,i1));
	dHj1  = squeeze(dH(:,j1,i1-1));
	dHj2  = squeeze(dH(:,j1,i1));
	Trj1(dHj1<1e-3)=0;
	Trj2(dHj2<1e-3)=0;
	U0(dHj1<1e-3) = 0;
	tr0   = 0.5*(Trj1+Trj2);
	dH0   = 0.5*(dHj1+dHj2);
	flx   = U0.*dH0.*tr0;
 	flx   = nx*U0.*dH0.*tr0*dd; % kg/m3*m*m/s*m=kg/s -> flux
      elseif dj==0
        Tri1  = squeeze(Tr(:,j1-1,i1));
	Tri2  = squeeze(Tr(:,j1,i1));
	V0    = squeeze(VN(:,j1,i1));
	dHi1  = squeeze(dH(:,j1-1,i1));
	dHi2  = squeeze(dH(:,j1,i1));
	Tri1(dHi1<1e-3)=0;
	Tri2(dHi2<1e-3)=0;
	V0(dHi1<1e-3)=0;
	tr0   = 0.5*(Tri1+Tri2);
	dH0   = 0.5*(dHi1+dHi2);
	flx   = ny*V0.*dH0.*tr0*dd; % kg/m3*m*m/s*m=kg/s -> flux
      end
      
      FLX(:,l) = flx;
      
    end;  % little segments along a section loop
%    keyboard
    
    TRFLX.Name = 'Atlantic OB';
    TRFLX.TM(nrec,1)=dnmb;
    TRFLX.FLX_kg_s(nrec,1)=nansum(nansum(FLX));

% Tracer mass within the region
    TRI = sub_intgr_tr(dH,Tr,Rmsk,Acell);
    TRFLX.SpongeReg_Mass_ton(nrec,1) = TRI.OverallMass_kg*1e-3; % mass inside sponge region

% Tracer mass over all domain
    TRI2 = sub_intgr_tr(dH,Tr,HH,Acell);
    TRFLX.Overall_Mass_ton(nrec,1) = TRI2.OverallMass_kg*1e-3; % mass inside whole domain

    dmm1=nansum(nansum(FLX));
    dmm2=TRI.OverallMass_kg*1e-3;
    dmm3=TRI2.OverallMass_kg*1e-3;
    rr  =dmm2/dmm3;
    fprintf('NetTrFlux %5.3d kg/s, Mass SpongeReg %7.5d ton, %6.3f%%\n',...
	    dmm1,dmm2,rr*100);
    fprintf('1 day processed %6.2f min \n\n',toc/60);
    
%    keyboard
    
    if s_mat>0 & mod(nrec,10)==0;
      fprintf('Saving %s\n',fmat);
      save(fmat,'TRFLX');
    end
    
  end % day
  if s_mat>0
    fprintf('Saving %s\n',fmat);
    save(fmat,'TRFLX');
  end
end

    
    

