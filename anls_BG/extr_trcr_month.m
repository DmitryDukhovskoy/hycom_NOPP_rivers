% This is an old code
% Bering Str coefficient is too hight
% depth-averaging is not accurate
%
% use extr_MassTrcr_month
%
% Save monthly mean tracer concentrations
% !!AVERAGED!! over some depth:
% 0 - 50 m - mixed layer
% 150 - 50 m
%
% river passive tracers 
% distributed along the Greenland Coast
%
% Adjust Bering Strait tracer concenctration:
% Distribute FW flux (0.08 Sv) over all points where Tracer=1
% # of j points = 7
% # of i points = 18
% total # of points (where HH>0 and above bottom)
% 123 points x 41 layers = 5043 grid points
% Thus tracer=1 unit is ~14-15.86 kg/m3 of tracer concentration
% 

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers
startup;

close all
clear


yr1=2015;
yr2=2015;


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
%cBr = 15.86; % coefficient for Bering Str. tracer
cBr = sub_BeringTrConc(0);% coefficient for Bering Str., tracer, kg/m3

%yr0 = 2005; 

lr1 = [0,-50];
lr2 = [-50,-150];
lr3 = [0,-10000]; % full water column

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
yrold = yr1;
dday = 5; 
for iyr = yr1:yr2
  yr=iyr;
%  for iday = 337:dday:366
%  id1=336;
  id1=1;
%  if yr == 2011,
%    id1 = datenum(2011,10,1)-datenum(2011,1,1)+1;
%  end
  
  for iday = id1:dday:365
    pthbin = sprintf('/nexsan/archive/%s_%3.3i/data/%i/',regn,expt,yr)

    dnmb=datenum(yr,1,1)+iday-1;
    DV=datevec(dnmb);
    imo=DV(2);
    fprintf('Processing Month %i\n',imo);

% skip processed months if needed
    if s_mat == 2
      fmat = sprintf('%strcr_dpthav_%i%2.2i.mat',pthmat,iyr,imo);
      if exist(fmat,'file')
	fprintf('Exist %s, skipping %i/%3.3i\n',fmat,yr,iday);
%    keyboard
	continue
      end
    end
    
    if mold~=imo
      if exist('TRCR','var') & s_mat>0 & mold>0
	nrec = TRCR(1).nrec;
	if nrec == 0,
	  error('# of saved records = 0');
	else
	  fprintf('Monthly mean, # of av. rcrds=%i\n',nrec);
	end
	
	fprintf('Saving TRCR...\n');
	for nTr=1:5
	  dmm = squeeze(TRCR(1).TR(nTr,:,:));
	  TRCR(1).TR(nTr,:,:) = dmm./nrec;
	  dmm = squeeze(TRCR(2).TR(nTr,:,:));
	  TRCR(2).TR(nTr,:,:) = dmm./nrec;
	  dmm = squeeze(TRCR(3).TR(nTr,:,:));
	  TRCR(3).TR(nTr,:,:) = dmm./nrec;
        end	
	
        fmat = sprintf('%strcr_dpthav_%i%2.2i.mat',pthmat,yrold,mold);
        fprintf('Saving %s\n',fmat);
        save(fmat,'TRCR');
      end
      
      clear TRCR
      TRCR(1).nrec=0;
      TRCR(2).nrec=0;
      TRCR(1).TR(5,mm,nn)=0;
      TRCR(2).TR(5,mm,nn)=0;
      TRCR(3).TR(5,mm,nn)=0;
      mold = imo;
      yrold = iyr;
      nrec = 0; 
    end

    fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
    finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);
  
    if ~exist(fina,'file');
      fprintf('Not found: %s\n\n',fina);
      continue;
    end
    

    nrec = nrec+1;
    tic;
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

% Convert m3/s -> kg/m3  
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
% Find z indices for specified layers:
%        iz1 = zeros(mm,nn);
%	iz2 = zeros(mm,nn);
%        for k=2:ll
%	  aa = squeeze(ZZb(k,:,:));
%	  I1 = find(aa>=-52);
%	  if ~isempty(I1),
%	    iz1(I1) = k;
%	  end
%	  I2 = find(aa>=-155);
%	  if ~isempty(I2)
%	    iz2(I2) = k;
%	  end
%	end;
%     end  % layer depths

% Integrate tracer over the depths:
% then divide by layer depths = mean tr. conc within 
% specified depth levels
      fprintf('Integrating Tracers over depths ...\n');
      dz1=zeros(mm,nn);
      dz2=zeros(mm,nn);
      dz3=zeros(mm,nn);
      dmm=zeros(mm,nn);
      smm1=zeros(mm,nn);
      smm2=zeros(mm,nn);
      smm3=zeros(mm,nn);
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
	
	I3 = find(~isnan(dZ));
	dz3(I3)=dz3(I3)+dZ(I3);
	smm3(I3) = smm3(I3)+squeeze(Ctr(k,I3))'.*dZ(I3);
%	fprintf('dz3=%6.1f\n',dz3(1500,1200));
	
      end      

      I1 = find(dz1>0);
      smm1(I1) = smm1(I1)./dz1(I1); % kg/m3
      
      I2 = find(dz2>0);
      smm2(I2) = smm2(I2)./dz2(I2);
      
      smm3 = smm3./dz3; % kg/m3
      
%      keyboard
%      dmm = squeeze(TRCR(1).TR(nTr,:,:));
      TRCR(1).TM = dnmb;
      TRCR(1).depth_av=lr1;
      TRCR(1).TR(nTr,:,:)=squeeze(TRCR(1).TR(nTr,:,:))+smm1;
      TRCR(1).nrec=nrec;
      TRCR(2).depth_av=lr2;
      TRCR(2).TR(nTr,:,:)=squeeze(TRCR(2).TR(nTr,:,:))+smm2;
      TRCR(2).nrec=nrec;
      TRCR(3).depth_av=lr3;
      TRCR(3).TR(nTr,:,:)=squeeze(TRCR(3).TR(nTr,:,:))+smm3;
      TRCR(3).nrec=nrec;
    end;   % nTr - tracers

%    if s_mat>0 
%      fmat = sprintf('%strcr_dpthav_%i%2.2i.mat',pthmat,iyr,imo);
%      fprintf('Saving %s\n',fmat);
%      save(fmat,'TRCR');
%    end
    fprintf('==== Processing 1 day, %i tracers %5.1f min\n',nTr,toc/60);
%    keyboard
  end % iday

% Save at the end of the year
  if exist('TRCR','var') & s_mat>0 
    nrec = TRCR(1).nrec;
    if nrec == 0,
      error('# of saved records = 0');
    else
      fprintf('Monthly mean, # of av. rcrds=%i\n',nrec);
    end

    fprintf('Saving TRCR...\n');
    for nTr=1:5
      dmm = squeeze(TRCR(1).TR(nTr,:,:));
      TRCR(1).TR(nTr,:,:) = dmm./nrec;
      dmm = squeeze(TRCR(2).TR(nTr,:,:));
      TRCR(2).TR(nTr,:,:) = dmm./nrec;
      dmm = squeeze(TRCR(3).TR(nTr,:,:));
      TRCR(3).TR(nTr,:,:) = dmm./nrec;
    end	

    fmat = sprintf('%strcr_dpthav_%i%2.2i.mat',pthmat,yrold,mold);
    fprintf('Saving %s\n',fmat);
    save(fmat,'TRCR');
  end

end

%keyboard



