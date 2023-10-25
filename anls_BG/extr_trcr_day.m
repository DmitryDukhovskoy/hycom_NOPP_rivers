% Save monthly mean tracer concentrations
% integrated over some depth:
% 0 - 50 m - mixed layer
% 150 - 50 m
%
% river passive tracers 
% distributed along the Greenland Coast
%
% Adjust Bering Strait tracer concenctration: <-- old, wrond adjustment
% see updated 04/2019
% However, this doesn't really matter 
% if the tracer is normalized by the overall tracer content 
% in the domain
% and then scaled by the FWFlux by the tracer source
%
% Distribute FW flux (0.08-0.1 Sv) over all points where Tracer=1
% 0.08 Sv = 2522.8 km3/yr
% 0.1 Sv = 3200 km3/yr
% Model long-term mean is 2778.9 km3/yr (anls_fluxes/calc_FWFlux_Bering.m)
%
%
% 

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers
startup;

close all
clear

s_mat = 1; % =0 - do not save mat file
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
%cBr = 14; % coefficient for Bering Str. tracer - not right
cBr = sub_BeringTrConc(0);% coefficient for Bering Str., tracer, kg/m3

%yr0 = 2005; 

lr1 = [0,-50];
lr2 = [-50,-150]; % note that in monthlym this is 0-150, also 1993-1998 daily

need to add overall depth for scaling
tracer concentration 
integrated over the whole domain volume


yr1=2016;
yr2=2016;

fprintf('Years to extract: %i-%i\n',yr1,yr2);

% Experiments:
% 102 - test simulation, HYCOM GLBb0.08 nest: 1 file kept constant (1993,1,1)

regn = 'ARCc0.08';
expt = 110;  
%pthfig  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/fig_trac/',expt);
pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);


ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

hmsk=HH;
hmsk(HH<0)=nan;

xlim1 = 0;
xlim2 = nn;
ylim1 = 0;
ylim2 = mm;

% Read fields:
ip1=1;
mold = 0;
dday = 7; 
for iyr = yr1:yr2
  TRCR = [];
  TRCR = struct;
  yr=iyr;
  fmat = sprintf('%strcr_dpthav_daily_%4.4i.mat',pthmat,iyr);
  id1=1;
  nrec = 0;
  for iday = id1:dday:365
    tic
    pthbin = sprintf('/nexsan/archive/%s_%3.3i/data/%i/',regn,expt,yr)

    dnmb=datenum(yr,1,1)+iday-1;
    DV=datevec(dnmb);
    imo=DV(2);
    fprintf('Processing Month %i\n',imo);

% skip processed days if needed
    if s_mat == 2
      if exist(fmat,'file')
	fprintf('Exist %s, skipping this year %i\n',fmat,yr);
	continue
      end
    end

    fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
    finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);
  
    if ~exist(fina,'file');
      fprintf('Not found: %s\n\n',fina);
      continue;
    end
    

    nrec = nrec+1;
    for nTr = 1:5
      fprintf('Tr# %i, %4.4i_%2.2i_%2.2i: %s\n',...
	  nTr,DV(1:3),fina);

      [F,nn,mm,ll] = read_hycom(fina,finb,'tracer','r_tracer',nTr);
      F(F>1e6)=nan;
      
% For Pacific Water,
% Convert nondim concentration 1 to 
% actual FW flux-based conc.
% used Woodgate & Aagaard estimate of Bering FW flux
%
      if nTr==5, 
	fprintf('Adjusting Bering Strait flux\n');
	F = F*cBr;
      end

% Average over upper 50 m
% Need layer thicknesses 
      Ctr = F; 

      if nTr == 1
	fprintf('Getting layer depths ...\n');
	fld='thknss';
	[F,n,m,l] = read_hycom(fina,finb,fld);
	F(F>1e10)=nan;
	F(F<0.1)=nan;
	F=F./rg;
	Lthck = F;
% Create Depth array of interface depths:
% Note these are BOTTOM interfaces 
% So Layer 1 is between interfaces 0m and ZZ(1)
        ZZb = F.*nan;
	ZZb(1,:,:) = 0;
	I = find(HH>=0);
	ZZb(1,I) = nan;
	for kk=1:ll
	  ZZb(kk+1,:,:)=ZZb(kk,:,:)-Lthck(kk,:,:);
	end
      end;
% Integrate tracer over the depths:
      fprintf('Integrating Tracers over depths ...\n');
      dz1=zeros(mm,nn);
      dz2=zeros(mm,nn);
      dmm=zeros(mm,nn);
      smm1=zeros(mm,nn);
      smm2=zeros(mm,nn);
      for k=1:ll
	aa = squeeze(ZZb(k+1,:,:));
	dZ = squeeze(abs(ZZb(k+1,:,:)-ZZb(k,:,:)));
	I1 = find(aa>=-52);
	if ~isempty(I1)
	  dmm = squeeze(Ctr(k,I1))'.*dZ(I1);
	  dz1(I1) = dz1(I1)+dZ(I1);
	  smm1(I1) = smm1(I1)+dmm;
	end
	
%	I2 = find(aa>=-155);
	I2 = find(aa>=-155 & aa<-52);
	if ~isempty(I2)
	  dmm = squeeze(Ctr(k,I2))'.*dZ(I2);
	  dz2(I2) = dz2(I2)+dZ(I2);
	  smm2(I2)= smm2(I2)+dmm;
	end
      end      

      I1 = find(dz1>0);
      smm1(I1) = smm1(I1)./dz1(I1); % kg/m3
      
      I2 = find(dz2>0);
      smm2(I2) = smm2(I2)./dz2(I2);
      
%      keyboard
%      dmm = squeeze(TRCR(1).TR(nTr,:,:));
      TRCR(nTr).TM(nrec) = dnmb;
      TRCR(nTr).depth_av1=lr1;
      TRCR(nTr).TR_lr1(nrec,:,:)=single(smm1);
      TRCR(nTr).depth_av2=lr2;
      TRCR(nTr).TR_lr2(nrec,:,:)=single(smm2);
    end;   % nTr - tracers

    fprintf('==== Processing 1 day %5.1f min\n',toc/60);
    
    if mod(nrec,10)==0 & exist('TRCR','var') & s_mat>0 
      fprintf('Saving TRCR...\n');

      fprintf('Saving %s\n',fmat);
      save(fmat,'TRCR','-v7.3');
    end
  
  end % iday
  
  if exist('TRCR','var') & s_mat>0 
    fprintf('Saving TRCR...\n');

    fprintf('Saving %s\n',fmat);
    save(fmat,'TRCR','-v7.3');
  end
end  % year

% Save at the end
if exist('TRCR','var') & s_mat>0 & mold>0
  if exist('TRCR','var') & s_mat>0 
    fprintf('Saving TRCR...\n');

    fmat = sprintf('%strcr_dpthav_daily_%4.4i.mat',pthmat,iyr);
    fprintf('Saving %s\n',fmat);
    save(fmat,'TRCR','-v7.3');
  end
end

%keyboard



