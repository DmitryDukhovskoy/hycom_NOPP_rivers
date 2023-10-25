% Save monthly mean tracer mass
% in every grid cell intergated over
% depth levels: 0-50, 50-150, 150-300, 0-btm
% Tracer is exactly integrated over the 
% specified depth levels, providing mass (kg)
%
% In order to get depth average Tr Conc within the 
% specified layers: Ctr(level=ilv)= Mass_Trcr(ilv)/(dZ*Area_cell)
% where dZ = Z1(ilv)-Z2(ilv) - layer thickness
%
% This is needed to assess S change 
% The codes is modified to
% extract & save 1 tracer at a time (nTr=1,...,5)
% Tracers: 1 - Greenland
%          2 - Mackenzie
%          3 - East Eurasian R.
%          4 - W. Eurasian R. 
%          5 - Bering Str. 
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers
startup;

close all
clear

nTr = 1; % tracer to extract
yr1=2003;
yr2=2004;
cmin = 1e-5; % treshold value of Tr. COnc.


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
fprintf('Bering Strait coeff=%5.2f\n',cBr);

%yr0 = 2005; 

% Vertical layers
LRS = load('LRS.dat');
nlrs= length(LRS); 

fprintf('Tracer: %i, Years to extract: %i-%i\n',nTr,yr1,yr2);

% Experiments:
% 102 - test simulation, HYCOM GLBb0.08 nest: 1 file kept constant (1993,1,1)

regn = 'ARCc0.08';
expt = 112;  
%pthfig  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/fig_trac/',expt);
pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
%pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
pthmat  = sprintf('/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/%3.3i/data_matTr/',expt);


ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

% Grid cell spacing
[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2

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
dday = 7; 
id1 = 1;
for iyr = yr1:yr2
  yr=iyr;
%  for iday = 337:dday:366
%  id1=336;
%  id1=337; % December
  
  for iday = id1:dday:365
    pthbin = sprintf('/nexsan/archive/%s_%3.3i/data/%i/',regn,expt,yr);

    if expt==112,
      pthbin=sprintf('/nexsan/hycom/ARCc0.08_%3.3i/data/%4.4i/',expt,yr);
    end

    dnmb=datenum(yr,1,1)+iday-1;
    DV=datevec(dnmb);
    imo=DV(2);
    fprintf('Processing Month %i\n',imo);

% skip processed months if needed
    if s_mat == 2
%      fmat = sprintf('%sMassTr_lrs_%i%2.2i.mat',pthmat,iyr,imo);
      fmat = sprintf('%sMassTr%2.2i_lrs_%i%2.2i.mat',pthmat,nTr,iyr,imo);
      if exist(fmat,'file')
	fprintf('Exist %s, skipping %i/%3.3i\n',fmat,yr,iday);
%    keyboard
	continue
      end
    end
    
% Month averaging    
    if mold~=imo
      if exist('TRCR','var') & s_mat>0 & mold>0
	nrec = TRCR(1).nrec;
	if nrec == 0,
	  error('# of saved records = 0');
	else
	  fprintf('Monthly mean, # of av. rcrds=%i\n',nrec);
	end
	
	fprintf('Saving TRCR...\n');
	for ik=1:nlrs
	  nrec = TRCR(ik).nrec;
	  if nrec==1, continue; end;
	  dmm = squeeze(TRCR(ik).MassTr_kg);
	  TRCR(ik).MassTr_kg = dmm./nrec;
          iii=sub2ind(size(HH),1500,1000);
	  trmm(ik)=TRCR(ik).MassTr_kg(iii);
	  TRCR(ik).nrec=1;
	  TRCR(ik).nrec_old=nrec;
	end
	
        trdm = nansum(nansum(TRCR(ik).MassTr_kg));
	for ik=1:nlrs-1
	  fprintf(' ===== Lr %i, pnt %i, Fraction rr=%8.6d\n',ik,iii,trmm(ik)/trdm);
	end

	dmm=TRCR(5).MassTr_kg(805:920,577:590);
	Im=find(isnan(dmm));
	if isnan(Im),
	  fprintf('NaN in Baffin Bay\n');
	  keyboard
	end

	fmat = sprintf('%sMassTr%2.2i_lrs_%i%2.2i.mat',pthmat,nTr,yrold,mold);
	if s_mat>0
          fprintf('Saving %s\n\n',fmat);
          save(fmat,'TRCR');
	end
%      keyboard
      end

      clear TRCR
      for ik=1:nlrs
	TRCR(ik).C_thrshold = cmin;
	TRCR(ik).Layer_Z1(1,1)=LRS(ik,1);
	TRCR(ik).Layer_Z2(1,1)=LRS(ik,2);
        TRCR(ik).nrec=0;
        TRCR(ik).MassTr_kg(mm,nn)=0;
%	TRCR(ik).layer_thck_m(mm,nn)=0;
      end
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
    tic;
%    for nTr = 1:5
    fprintf('nrec= %i, Tr# %i, %4.4i_%2.2i_%2.2i: %s\n',...
	  nrec, nTr,DV(1:3),fina);
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

    Ctr = F; % kg/m3
    Ctr(Ctr<cmin)=nan;

%      if nTr == 1
    fprintf('Getting layer depths ...\n');
    [ZM,ZZ] = sub_zz_zm(fina, finb,HH,'f_btm',1);
%      end
% Integrate tracer over the depths:
% to get tracer mass (kg)
% specified depth levels
    fprintf('Integrating Tracers over depths ...\n');
% Integrate S over layer ilv
    for ilv=1:nlrs
      zz1=LRS(ilv,1);
      zz2=LRS(ilv,2);
      dz=abs(zz2-zz1);
      zbtm=-9000;
      MTr = sub_zLayer_average(HH,ZZ,Ctr,zbtm,zz1,zz2); % note returned is tr conc kg/m3!
      if isempty(MTr), error('Sav is empty'), end;
      if zz2>zbtm % av. in layer
	hZ=HH*0+dz;
	hZ(HH>=0)=nan;
      else
	hZ=abs(HH);
	hZ(HH>=0)=nan;
      end

      MTr = MTr.*hZ.*Acell; % Layer-averaged TrConc -> Mass, kg
%      mtrB = MTr;
%      mtr  = MTr;
%      mtrB0= MTr; % tot mass without any Ctr thresholding
%keyboard      
      MTr(MTr<=0)=nan;

      tmm = squeeze(TRCR(ilv).MassTr_kg);
      TRCR(ilv).TM = dnmb;
      TRCR(ilv).depths = LRS(ilv,:);
      TRCR(ilv).MassTr_kg = tmm+MTr;
      TRCR(ilv).nrec=nrec;

    end % levels

%    end;   % nTr - tracers

    fprintf('==== Processing 1 day, %i tracers %5.1f min\n',nTr,toc/60);

  end % iday

% Save at the end of last year
  if exist('TRCR','var') & s_mat>0 & iyr==yr2
    nrec = TRCR(1).nrec;
    if nrec == 0,
      error('# of saved records = 0');
    else
      fprintf('Monthly mean, # of av. rcrds=%i\n',nrec);
    end

    fprintf('End of Last Year in the Loop: Saving TRCR...\n');
%    for nTr=1:5
      for ik=1:nlrs
	nrec = TRCR(ik).nrec;
	if nrec==1, continue; end;
	dmm = squeeze(TRCR(ik).MassTr_kg);
	TRCR(ik).MassTr_kg = dmm./nrec;
	i=2265680;
	trmm(ik)=TRCR(ik).MassTr_kg(i);
	TRCR(ik).nrec=1;
	TRCR(ik).nrec_old=nrec;
      end
      trdm = nansum(nansum(TRCR(ik).MassTr_kg));
      for ik=1:nlrs-1
	fprintf(' ===== Lr %i, pnt %i, Fraction rr=%8.6d\n',ik,i,trmm(ik)/trdm);
      end

      
%    fmat = sprintf('%sMassTr_lrs_%i%2.2i.mat',pthmat,yrold,mold);
    fmat = sprintf('%sMassTr%2.2i_lrs_%i%2.2i.mat',pthmat,nTr,yrold,mold);
    fprintf('  ====   Saving %s\n\n',fmat);
%    keyboard
    save(fmat,'TRCR');
  end

end

%keyboard



