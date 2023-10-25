% Compute volume integrated mass of the tracer
% for specified regions
% save for whole water column
%
% VTRCR - structured array with 
% vertical profiles of tracer mass
% in HYCOM Layers averaged over each region
% for nbx regions
% 
%
% monthly mean
% within the specified domain IN
% Use only deep basin - deeper than 500m
% If IN is empty - the whole domain
% nTr - tracer # that is integrated
% Acell - grid cell area, m2
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

regn = 'ARCc0.08';
expt = 112;  

nvl   = 41;  % # of v. layers in the model
nTr   = 1;   % tracer to plot
s_mat = 1; % =1 - save mat file flag
           % =2 - load saved fields
% 1993        1999        2005     2011
YR1 = 2008;
YR2 = 2009;
dyr = 1;
h0  = -800; % cutoff depth, m

rg = 9806;

fprintf('Tracer: %i, Depths<%4.1fm, Layer integrated, %i:%i:%i\n',...
	nTr,h0,YR1,dyr,YR2);

% Experiments:
pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc0.08/%3.3i/data_mat/',expt);


%ftopo = sprintf('%s/depth_%s_17DD.nc',pthtopo,regn); % 
ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(HH);
[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2
MVOL = [];
Hmsk = HH;
Hmsk(HH>h0)=100;


% Define Regions:
f_pltbox=0;
BX = sub_deep_regions(HH,LON,LAT,h0,f_pltbox);
%BX = sub_define_boxes(HH,LON,LAT,f_pltbox);
%nbx = length(BX);
nbx = length(BX); % 

if s_mat==1
  fprintf('Mat file will be saved\n');
else
  fprintf('Mat file WILL NOT be saved\n');
end  

mold  = 0;
yrold = YR1;
dday  = 14; 
for iyr=YR1:dyr:YR2
  yr   = iyr;
%  fmat = sprintf('%s%3.3i_VrtLrMass_mo_Tr%2.2i_%i.mat',...
%		 pthmat,expt,nTr,iyr);
  fmat = sprintf('%strcr%2.2i_regn_VrtLrMass_%i.mat',pthmat,nTr,iyr);
    
  cc=0;
  for ibx=1:nbx
    VTRCR(ibx).Title= sprintf('%s Tracer Mass integrated over boxes by layers',regn);
    VTRCR(ibx).info = sprintf(' Inside isobath %5.1f',h0);
    VTRCR(ibx).Name = BX(ibx).Name;
    VTRCR(ibx).XY   = BX(ibx).XY;
    VTRCR(ibx).IJ   = BX(ibx).IJ;
    VTRCR(ibx).IN   = BX(ibx).IN;
  end
  
  for iday=1:dday:365
    dnmb=datenum(yr,1,1)+iday-1;
    DV=datevec(dnmb);
    imo=DV(2);
    
    if mold~=imo
      if exist('VTRCR','var') & s_mat>0 & mold>0
	nrec = VTRCR(ibx).nrec;
	if nrec == 0,
	  error('# of saved records = 0');
	else
	  fprintf('Monthly mean, # of av. rcrds=%i\n',nrec);
	end
	
	fprintf('Saving VTRCR...\n');
	for ibx=1:nbx
	  dmm = squeeze(VTRCR(ibx).TR(:,mold));
	  dzm = VTRCR(ibx).DZM(:,mold);
%	  VTRCR(ibx).TR(:,mold)  = dmm./dzm;  % average tracer mass in layers, kg
	  VTRCR(ibx).TR(:,mold)  = dmm./nrec;  % average tracer mass in layers, kg
	  VTRCR(ibx).DZM(:,mold) = dzm./nrec; % average layer thkn, m
% check layer thicknesses - should be comparable to total depth	  
	  gmm = dzm./nrec;
	  hdp = sum(gmm);
	  mh0  = VTRCR(ibx).meanH;
	  epsl = abs(1-abs(mh0/hdp));
	  
	  fprintf('Estimated total depth=%6.1fm, mean depth=%6.1f\n',...
		  -hdp, mh0);
	  if epsl>0.02,
	    fprintf('Check dzm layer thicknesses, do not match total H\n');
	    keyboard;
	  end
	  
        end	

	
        fprintf('Saving %s\n',fmat);
        save(fmat,'VTRCR');
%	keyboard
      end
      
      for ibx=1:nbx 
	VTRCR(ibx).nrec         = 0;
	VTRCR(ibx).TR (1:nvl,imo) = 0;
	VTRCR(ibx).DZM(1:nvl,imo) = 0;
      end
      mold  = imo;
      yrold = iyr;
      nrec  = 0; 
    end
   
%pthbin = sprintf('/nexsan/hycom/ARCc0.08_011/data/%i/',yr);  
    pthbin = sprintf('/nexsan/archive/%s_%3.3i/data/%i/',regn,expt,yr);
    if expt==112
      pthbin = sprintf('/nexsan/hycom/ARCc0.08_112/data/%i/',yr);
    end
    fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
    finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);

    fprintf('\n:::: Analyzing %i/%2.2i/%2.2i   ::::\n\n',DV(1:3));
    cc = cc+1;
    nrec = nrec+1;
    
    tic;
  
%    fprintf('Computing vol-integrated Tracer mass ...\n');
    [F,n,m,l] = read_hycom(fina,finb,'tracer','r_tracer',nTr);
    F(F>1e6)=nan;
    Tr=F;

    [F,n,m,l] = read_hycom(fina,finb,'thknss');
    F=F./rg;
    F(F>1e10)=0;
    F(F<1e-2)=0;
    dH=F; 
    
    for ibx = 1:nbx
      IN = BX(ibx).IN;
      Hmsk=HH*nan;
      Hmsk(IN)=HH(IN);
      tri = sub_intgr_tr(dH,Tr,Hmsk,Acell);
      dmm = VTRCR(ibx).TR(:,imo);
%      VTRCR(ibx).TR(:,imo) = dmm+tri.Vert_LrMass_kg.*tri.Mean_LThkn;
      VTRCR(ibx).TR(:,imo) = dmm+tri.Vert_LrMass_kg; % mass in each layer
      zmm = VTRCR(ibx).DZM(:,imo);
      VTRCR(ibx).DZM(:,imo) = zmm+tri.Mean_LThkn;
      VTRCR(ibx).nrec = nrec;
      VTRCR(ibx).TM(nrec,imo) = dnmb;
      VTRCR(ibx).meanH = nanmean(HH(IN));
    end
    
    fprintf('1 day processed %6.2fmin\n\n',toc/60);
%    keyboard
    
  end
  
  if s_mat==1 
    for ibx=1:nbx
      nrec = VTRCR(ibx).nrec;
      if nrec == 0, continue; end; 
      dmm = squeeze(VTRCR(ibx).TR(:,mold));
      dzm = VTRCR(ibx).DZM(:,mold);
%      VTRCR(ibx).TR(:,mold)  = dmm./dzm;  % average tracer mass in layers, kg
      VTRCR(ibx).TR(:,mold)  = dmm./nrec;  % average tracer mass in layers, kg
      VTRCR(ibx).DZM(:,mold) = dzm./nrec; % average layer thkn, m
% check layer thicknesses - should be comparable to total depth	  
      gmm = dzm./nrec;
      hdp = sum(gmm);
      mh0  = VTRCR(ibx).meanH;
      epsl = abs(1-abs(mh0/hdp));

      fprintf('== Estimated total depth=%6.1fm, mean depth=%6.1f\n',...
	      -hdp, mh0);
      if epsl>0.02,
	fprintf('== Check dzm layer thicknesses, do not match total H\n');
	keyboard;
      end
      
      VTRCR(ibx).nrec = 0;
	  
    end	
    fprintf('========  End of year: Saving %s\n',fmat);
    save(fmat,'VTRCR');
    
    mold = 0; % to prevent saving it 2nd time
    clear VTRCR
  end

end; % time loop

%exit
