% 0.04 HYCOM2.3/2.3-CICE5 GOFS3.5 
% 
%
% Save monthly mean U,V fields
% integrated over some depth (tansports m3/s):
% 0 - 50 m - mixed layer
% Use Gauss theorem to estimate divergnece
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

yr1=2017;
yr2=yr1;  % annual fields are saved


zz0 = 50;
iz0 = [];
dgr = 31; % # of grid cells (dx=dy=dgr) to calculate 
         % the contour integral of normal U
	 % use odd number 1,3,5,...
dij = (dgr-1)/2; % (i,j) +/-dij

s_mat = 2; % =0 - do not save mat file
           % =1 - save mat
      	   % =2 - skip months where mat file exist

if s_mat==0,
  fprintf('Mat file is not created\n');
elseif s_mat == 1
  fprintf('Mat file will be saved, old mat file will be overridden\n');
elseif s_mat == 2
  fprintf('Extraction is skipped for months where old mat files exist\n');
end


rg = 9806; 
lr1 = [0,-zz0];

fprintf('Years to extract: %i-%i\n',yr1,yr2);

% Choose experiment:
ixx    = 9; % experiment name and dir - check with EXPT - expt 023
EXPT   = sub_cice_experiments;
expt   = EXPT(ixx).Nmb;
texpt  = EXPT(ixx).cice_opt; % CICE options for sens. experiments
res    = EXPT(ixx).res;

regn = 'ARCc0.04';
rg = 9806;
hgg=1e20;


pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';
pthmat = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/%3.3i/Uvort/',expt);

ftopo = sprintf('%s/depth_%s_17DD.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);
Acell = DX.*DY;

Ideep = find(HH<-500);

% Spatial filtering:
pgrd = 101;  % must be odd
Hmsk = HH;
Hmsk(HH<0)=1;
Hmsk(HH>=0)=0;
Hmsk(3840:end,:)=0;
Hmsk(200:1000,2500:end)=0;
Hmsk(900:1300,2560:2760)=0;
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
  id1=2;
  
  fmat = sprintf('%s%3.3i_divU%3.3im_%i.mat',...
		   pthmat,expt,abs(zz0),iyr);

   IDD = [id1:dday:365];
   TM0 = datenum(yr,1,1)+IDD-1;
   DV0 = datevec(TM0);
   dnmb0=0;
% For saved fields, find last record
  if s_mat==2
    fprintf('Loading saved %s\n',fmat);
    load(fmat);
 
    mlast = length(DIVU); 
    ilast = max(find(DV0(:,2)==mlast));
    dnmb0 = TM0(ilast);
   
    fprintf('Last saved field: %s\n',datestr(dnmb0));     
 
  end

  for iday = IDD
    dnmb=datenum(yr,1,1)+iday-1;
    if dnmb<=dnmb0; 
      fprintf('===>  Field allready processed, skipping %s ...\n',datestr(dnmb))
      continue; 
    end

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

% Read fields:
    pthbin = sprintf('/nexsan/people/ddmitry/hycom/ARCc0.04_%3.3i/data/%4.4i_%s/',...
                    expt,yr,texpt);
% For expt 0.23 use first 3 months from 02.2 - 02.3 started from 02.2
    expt   = EXPT(ixx).Nmb;
    if ixx == 9 & dnmb<datenum(2017,4,1),
      pthbin = '/nexsan/people/ddmitry/hycom/ARCc0.04_022/data/2017_BL99Tfrz/';
      expt = 22;
    end

%    fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
%    finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);
    fina = sprintf('%s%3.3i_archv.%4.4i_%3.3i_00.a',pthbin,expt,yr,iday);
    finb = sprintf('%s%3.3i_archv.%4.4i_%3.3i_00.b',pthbin,expt,yr,iday);
  
    if ~exist(fina,'file');
      fprintf('Not found: %s\n\n',fina);
      continue;
    end
    

    nrec = nrec+1;
    DIVU(mold).nrec=nrec;
    tic;
    
    [F,n,m,nlr] = read_hycom(fina,finb,'u-vel.');
    F(F>1e6)=nan;
    U=squeeze(F);

		[F,n,m,llr] = read_hycom(fina,finb,'u_btrop');
		F(F>hgg)=0;
		Ub=squeeze(F);
%
% For archv velocities need to add barotropic comp.
		for ilr=1:nlr
			U(ilr,:,:)=squeeze(U(ilr,:,:))+Ub;
		end

    [F,n,m,l] = read_hycom(fina,finb,'v-vel.');
    F(F>1e6)=nan;
    V=squeeze(F);

		[F,n,m,llr] = read_hycom(fina,finb,'v_btrop');
		F(F>hgg)=0;
		Vb=squeeze(F);

		for ilr=1:nlr
			V(ilr,:,:)=squeeze(V(ilr,:,:))+Vb;
		end

    fprintf('Getting layer depths ...\n');
    fld='thknss';
    [F,n,m,ll] = read_hycom(fina,finb,fld);
    F(F>1e10)=nan;
    F(F<0.1)=nan;
    F=F./rg;
    dH = F;

%    fprintf('Getting layer depths ...\n');
%    [ZM,ZZ] = sub_zz_zm(fina, finb,HH,'f_btm',1);

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
				dh1=squeeze(dH(k,:,ii-1));
				dh2=squeeze(dH(k,:,ii));
				dh1(isnan(dh1))=0;
				dh2(isnan(dh2))=0;
				dHu(k,:,ii)=0.5*(dh1+dh2);
      end

      for jj=2:mm
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
    divU = sub_gauss_divU(Hmsk,Utr,Vtr,dij);
    sdivU = sdivU+divU;
    
%    
% ==============================    
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!1    
    f_chck = 0;
    if f_chck==1
      i=646;
      j=491;
      
% Note that ut and vt are not collocated
      aa = Utr(j-dij:j+dij,i-dij:i+dij);
      bb = Vtr(j-dij:j+dij,i-dij:i+dij);
      figure(10); clf;
      quiver(aa,bb);
      hold
      plot(dij+1,dij+1,'r*');

      
      ut=Utr(j-dij:j+dij,i-dij:i+dij+1);
      vt=Vtr(j-dij:j+dij+1,i-dij:i+dij);
      ut(isnan(ut))=0;
      vt(isnan(vt))=0;
      iUn = -sum(ut(:,1))+sum(ut(:,end))-sum(vt(:,1))+sum(vt(:,end));
      
% Gauss theorem for calculating mean divergence:
% This doesn't seem to work in HYCOM
% archm fields - lateral integrated flux
% does not equal volume integrated divU
%
% Gauss theorem: surface integral - lateral boundary:
      Und = zeros(iz0,1);
      for ka=1:iz0
				fprintf('ka=%i \n',ka);
				ut  = squeeze(U(ka,:,:));
%        ut  = squeeze(U(ka,j-dij:j+dij,i-dij:i+dij+1));
        ut  = squeeze(uu(j-dij:j+dij,i-dij:i+dij+1));
				dhu = squeeze(dHu(ka,:,:));
				dhu = squeeze(dhu(j-dij:j+dij,i-dij:i+dij+1));
				dyy = DY(j-dij:j+dij,i-dij:i+dij+1);
	
        vt  = squeeze(V(ka,:,:));
        vt  = squeeze(vt(j-dij:j+dij+1,i-dij:i+dij));
				dhv = squeeze(dHv(ka,:,:));
				dhv = squeeze(dhv(j-dij:j+dij+1,i-dij:i+dij));
				dxx = DX(j-dij:j+dij+1,i-dij:i+dij);

				ut(isnan(ut))=0;
				vt(isnan(vt))=0;
				dhu(isnan(dhu))=0;
				dhv(isnan(dhv))=0;

				dmm = -sum(ut(:,1).*dhu(:,1).*dyy(:,1))+...
							 sum(ut(:,end).*dhu(:,end).*dyy(:,end))-...
							 sum(vt(:,1).*dhv(:,1).*dxx(:,1))+...
							 sum(vt(:,end).*dhv(:,end).*dxx(:,end));
				Und(ka,1) = Und(ka)+dmm;
      end
      sUnd = sum(Und); % should equal iUn and ideally sdvU	
      
	
% volume integral of divergence;
% This doesn't seem to work well
      dvU = zeros(iz0,1);
      hha = 0;
      for ka=1:iz0
				ua = squeeze(U(ka,:,:));
				va = squeeze(V(ka,:,:));
				dha= squeeze(dH(ka,:,:));
				hha = hha+dha(j,i);
				ua(isnan(ua))=0;
				va(isnan(va))=0;
				dha(isnan(dha))=0;
				fprintf('ka=%i, z=%6.1fm\n',ka,-hha);
        for ia=i-dij:i+dij
					for ja=j-dij:j+dij
						dx = DX(ja,ia);
						dy = DY(ja,ia);
						divu = (ua(ja,ia+1)-ua(ja,ia))/dx+(va(ja+1,ia)-va(ja,ia))/dy;
						dvU(ka,1)  = dvU(ka)+divu*dx*dy*dha(ja,ia);
					end
        end
      end
      sdvU = sum(dvU);
    end
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!1    
% ==============================    
    
    fprintf('==== Processing 1 day, %5.1f min\n',toc/60);
%    keyboard
  end % iday


end

if s_mat>0 & mold>0
  nrec = DIVU(mold).nrec;
  if nrec == 0, return; end;


  divU = sdivU/nrec;
  
  DIVU(mold).nrec = nrec;
  DIVU(mold).NgrdPnts_dx = dgr;
  DIVU(mold).divU_m3_sec = divU;
  DIVU(mold).SurfLayer_integrated_m  = zz0;

%        fmat = sprintf('%s%3.3i_UVintgr%3.3im_%i.mat',...
%		     pthmat,expt,abs(zz0),iyr);

  fprintf('===== END of YEAR ====  Saving %s\n',fmat);
  save(fmat,'DIVU');
  
  nrec  = 0; 
  DIVU(mold).nrec = 0;
end

  

  



