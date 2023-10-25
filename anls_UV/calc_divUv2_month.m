% Save monthly mean U,V fields
% integrated over some depth (tansports m3/s):
% 0 - 50 m - mixed layer
% Approximate divU = du/dx+dv/dy
% as vol fluxes across grid sides: dlt(Tx)+dlt(Ty)
% Tx = intgr(z0:0) u*dh*dy - m3/s
% Thus this is Gauss theorem but for 1 grid, gives
% m3/s - volume flux divergence, i.e. how much
% water is pumped down (if divU<0) or succed up (divU>0)
% 
%            V(i,j+1)
%   --------|------
%   |             |
%   |             |
% U -      *(i,j) - U(i+1,j)
%   |             |
%   |             |
%   --------|------
%          V(i,j)
%
% Do not collocate U and V!
% But need to interpolate dH into U,V grid points
%

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers
startup;

close all
clear

yr1=2013;
yr2=2013;

zz0 = 50;
iz0 = [];
dgr = 15; % # of grid cells (dx=dy=dgr) to calculate 
         % the contour integral of normal U
	 % use odd number 1,3,5,...
dij = (dgr-1)/2; % (i,j) +/-dij

s_mat = 2; % =0 - do not save mat file
           % =1 - save tracer and overwrite existing mat
	   % =2 - skip months where mat file exist

if s_mat==0,
  fprintf('Mat file is not created\n');
elseif s_mat == 1
  fprintf('Mat file will be saved, old mat file will be overridden\n');
elseif s_mat == 2
  fprintf('Extraction is skipped for months where old mat files exist\n');
end


rg = 9806; 
%cBr = 15.86; % coefficient for Bering Str. tracer
lr1 = [0,-zz0];

fprintf('Years to extract: %i-%i\n',yr1,yr2);

% Experiments:
% 102 - test simulation, HYCOM GLBb0.08 nest: 1 file kept constant (1993,1,1)

regn = 'ARCc0.08';
expt = 110;  
%pthfig  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/fig_trac/',expt);
pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_mat2/';

ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);
Acell = DX.*DY;

Ideep = find(HH<-500);

% Spatial filtering:
pgrd = 51;
Hmsk = HH;
Hmsk(HH<0)=1;
Hmsk(HH>=0)=0;
Hmsk(1920:end,:)=0;
Hmsk(100:500,1250:end)=0;
Hmsk(450:650,1280:1380)=0;
Hmsk(1:pgrd+1,:) = 0;
Hmsk(mm-pgrd:mm,:) = 0;
Hmsk(:,1:pgrd+1) = 0;
Hmsk(:,nn-pgrd:nn) = 0;
Imsk = find(Hmsk==1);

xlim1 = 0;
xlim2 = nn;
ylim1 = 0;
ylim2 = mm;

% Read fields:
ip1=1;
mold = 0;
yrold = yr1;
dday = 6; 
for iyr = yr1:yr2
  yr=iyr;
  id1=1;
  
  fmat = sprintf('%s%3.3i_divU%3.3im_v2_%i.mat',...
		   pthmat,expt,abs(zz0),iyr);

  if s_mat==2
    if ~exist(fmat,'file');
      id1=1;
    else
      fprintf('::: Found existing file %s\n',fmat);
      load(fmat);
      nmm = length(DIVU);
      fprintf('::: %i month done\n',nmm);
      if nmm==12, continue; end;
      mold=0;
      dmn = datenum(iyr,nmm+1,1);
      ndays = dmn-datenum(iyr,1,1)+1;
      nrr = floor((ndays-1)/dday)+1;
      id1 = 1+nrr*dday;
      dnew = datenum(iyr,1,1)+id1;
      fprintf('::: Starting from day %s, id1=%i\n\n',datestr(dnew),id1);
    end
%    keyboard
  end
  
  for iday = id1:dday:365
    pthbin = sprintf('/nexsan/archive/%s_%3.3i/data/%i/',regn,expt,yr)

    dnmb=datenum(yr,1,1)+iday-1;
    DV=datevec(dnmb);
    imo=DV(2);
    fprintf('Processing %i/%2.2i/%2.2i\n',DV(1:3));

% ================== 
% When current month is done:
% Average div fields
% Filter (spatial filtering
% with 2D Gaussian)
% ==================    
    if mold~=imo
      if exist('DIVU','var') & s_mat>0 & mold>0
	nrec = DIVU(mold).nrec;
	if nrec == 0,
	  error('# of saved records = 0');
	else
	  fprintf('Monthly mean, # of av. rcrds=%i\n',nrec);
	end
	
%	Utr = Utr./nrec; % U transp intgr in upper zz0 m
%	Vtr = Vtr./nrec; % V transp intgr in upper zz0 m
	
% Sum transports:
% over Y faces:
% Note orientation of the normal unit vector
% directed outside the contour
%        divU = sub_gauss_divU(Hmsk,Utr,Vtr,dij);

        divU = sdivU/nrec;
        divU = sub_fltr(divU,pgrd,Hmsk);
	
	DIVU(mold).nrec = nrec;
	DIVU(mold).Info1 = 'Divergence, m3/s, mean from Gauss Theorem';
	DIVU(mold).Info2 = sprintf('Monthly mean, Filtered %i pnts, Gaussian',pgrd);
	DIVU(mold).NgrdPnts_dx = dgr;
	DIVU(mold).divU_m3_sec = divU;
	DIVU(mold).SurfLayer_integrated_m  = zz0;
	
        fprintf('Saving %s\n',fmat);
        save(fmat,'DIVU');
	
%	keyboard
      end
      
      Utr = HH*0;
      Vtr = HH*0;
      sdivU = HH*0;
      
      DIVU(imo).nrec    = 0;
      DIVU(imo).divU_m3_sec(mm,nn)= 0;
      mold  = imo;
      yrold = iyr;
      nrec  = 0; 
    end

    fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
    finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);
  
    if ~exist(fina,'file');
      fprintf('Not found: %s\n\n',fina);
      continue;
    end
    

    nrec = nrec+1;
    DIVU(mold).nrec=nrec;
    tic;
    
    [F,n,m,l] = read_hycom(fina,finb,'u-vel.');
    F(F>1e6)=nan;
    U=squeeze(F);

    [F,n,m,l] = read_hycom(fina,finb,'v-vel.');
    F(F>1e6)=nan;
    V=squeeze(F);

    fprintf('Getting layer depths ...\n');
    fld='thknss';
    [F,n,m,ll] = read_hycom(fina,finb,fld);
    F(F>1e10)=nan;
    F(F<0.1)=nan;
    F=F./rg;
    dH = F;

% Define layer that is ~=zz0  
% assume little seasonal and interannual
% variabiltiy at layer depths at 50m 
% - mostly geopotential layers
% thus do not need to find iz0 for every record
    if (isempty(iz0))
% Create Depth array of interface depths:
      ZZb = F.*nan;
      ZZb(1,:,:) = 0;
      I = find(HH>=0);
      ZZb(1,I) = nan;
      for kk=1:ll
	ZZb(kk+1,:,:)=ZZb(kk,:,:)-dH(kk,:,:);
      end

      for k=1:ll
	aa = squeeze(ZZb(k,:,:));
	dmn= nanmean(aa(Ideep));
	iz0 = k;
	if abs(dmn)>abs(zz0); break; end;
      end
      fprintf(' Integrating down to %i Layer\n',iz0);
    end
    
 %   Tx=U*0;
 %   Ty=V*0;
    fprintf('dHu, dHv, ...\n');
    dHu = U*0; % collocated in U points
    dHv = V*0;
    for k=1:iz0
      for ii=2:nn
   %     u1 = U(:,ii);
   %     u2 = U(:,ii+1);
   %     uu = 0.5*(u1+u2);
	dh1=squeeze(dH(k,:,ii-1));
	dh2=squeeze(dH(k,:,ii));
	dh1(isnan(dh1))=0;
	dh2(isnan(dh2))=0;
	dHu(k,:,ii)=0.5*(dh1+dh2);
    %    Tx(:,ii)=uu.*dh;
  %      Tx(:,ii) = uu; % collocated U in the layer
      end

      for jj=2:mm
  %      v1 = V(jj,:);
  %      v2 = V(jj+1,:);
  %      vv = 0.5*(v1+v2);
  %      Ty(jj,:)=vv;
	dh1=squeeze(dH(k,jj-1,:));
	dh2=squeeze(dH(k,jj,:));
	dh1(isnan(dh1))=0;
	dh2(isnan(dh2))=0;
	dHv(k,jj,:)=0.5*(dh1+dh2);
      end
    end
    
  
% Calculate transports (not collocated!)
    Utr = HH*0;
    Vtr = HH*0;
    for k=1:iz0
      dhu = squeeze(dHu(k,:,:));
      dhv = squeeze(dHv(k,:,:));
      uu  = squeeze(U(k,:,:));
      vv  = squeeze(V(k,:,:));
      utr = uu.*dhu.*DY; % transport, m3/s
      vtr = vv.*dhv.*DX; % transport, m3/s
      Utr = Utr+utr; % integrate over iz0 layers and month
      Vtr = Vtr+vtr; % integrate over iz0 layers and month
%      DHU = DHU+dhu;
%      DHV = DHV+dhv;
    end
      
    dUdX = HH*0;
    dVdY = HH*0;
    for ii=1:nn-1
      uL = Utr(:,ii);
      uR = Utr(:,ii+1);
      dx = DX(:,ii);
      dUdX(:,ii) = (uR-uL); % m3/s
    end
    for jj=1:mm-1
      vU = Vtr(jj+1,:);
      vB = Vtr(jj,:);
      dy = DY(jj,:);
      dVdY(jj,:) = (vU-vB);
    end

    divU = dUdX+dVdY; % m3/s - integrated mass flux over lateral faces
    sdivU = sdivU+divU;
    
    fprintf('==== Processing 1 day, %5.1f min\n',toc/60);
%    keyboard
  end % iday


end

if s_mat>0 & mold>0
  nrec = DIVU(mold).nrec;
  if nrec == 0, return; end;


  divU = sdivU/nrec;
  divU = sub_fltr(divU,pgrd,Hmsk);
  
  DIVU(mold).nrec = nrec;
  DIVU(mold).NgrdPnts_dx = dgr;
  DIVU(mold).divU_m3_sec = divU; % 
  DIVU(mold).SurfLayer_integrated_m  = zz0;

%        fmat = sprintf('%s%3.3i_UVintgr%3.3im_%i.mat',...
%		     pthmat,expt,abs(zz0),iyr);

  fprintf('===== END of YEAR ====  Saving %s\n',fmat);
  save(fmat,'DIVU');
  
  nrec  = 0; 
  DIVU(mold).nrec = 0;
end

  

  



