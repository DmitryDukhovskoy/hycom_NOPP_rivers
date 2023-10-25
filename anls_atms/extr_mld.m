% Extract MLD 
% Output from the model
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

yr1=1997;
yr2=1997;

fprintf('Years to extract: %i-%i\n',yr1,yr2);

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

regn = 'ARCc0.08';
expt = 110; 
rg   = 9806; 
%pthfig  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/fig_trac/',expt);
pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);

ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

btx = 'surf_heat_flx.m';


% Read fields:
ip1=1;
mold = 0;
yrold = yr1;
dday = 5; 
for iyr = yr1:yr2
  yr=iyr;
  id1=1;
  
  for iday = id1:dday:365
    pthbin = sprintf('/nexsan/archive/%s_%3.3i/data/%i/',regn,expt,yr);

    dnmb=datenum(yr,1,1)+iday-1;
    DV=datevec(dnmb);
    imo=DV(2);

% skip processed months if needed
    if s_mat == 2
      fmat = sprintf('%smld_%i.mat',pthmat,yrold);
      if exist(fmat,'file')
	fprintf('==> Exist %s skipping %i/%3.3i\n',fmat,yr,iday);
	continue
      end
    end
    
    tic;
    if mold~=imo
      if exist('MLD','var') & s_mat>0 & mold>0
	nrec = MLD(mold).nrec;
	if nrec == 0,
	  error('# of saved records = 0');
	else
          fprintf('Mo=%i, Monthly mean, # of av. rcrds=%i\n',mold,nrec);
	end
	
	dmm = squeeze(MLD(mold).mld);
	MLD(mold).mld = dmm./nrec;

	if mold==12
          fmat = sprintf('%smld_%i.mat',pthmat,yrold);
          fprintf('Saving %s\n',fmat);
          save(fmat,'MLD');
	end
	
      end
      if yrold~=iyr
        clear MLD
      end
      MLD(imo).nrec=0;
      MLD(imo).mld(mm,nn)=0;
      mold = imo;
      yrold = iyr;
      nrec = 0; 
      
    end

    fprintf('Processing Month %i\n',imo);
    fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
    finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);
  
    if ~exist(fina,'file');
      fprintf('!!!!!   Not found: %s  !!!!\n\n',fina);
      continue;
    end
    

    nrec = nrec+1;
    fprintf('Reading %4.4i_%2.2i_%2.2i: %s\n',DV(1:3),fina);
    [F,nn,mm,ll] = read_hycom(fina,finb,'mix_dpth');
    F=-1*(squeeze(F)./rg);
    F(F>1e12)=nan;
    F(HH>0)=nan;
    
    MLD(imo).TM   = dnmb;
    MLD(imo).nrec = nrec;
    MLD(imo).mld = MLD(imo).mld+F;
    
    fprintf('==== Process 1 day: %6.4f min\n\n',toc/60);
%   keyboard
  end % day

end % year

if exist('MLD','var') & s_mat>0
  nrec = MLD(mold).nrec;
  if nrec == 0,
    fprintf('# of saved records = 0, mat file not saved');
    return
  else
    fprintf('Mo=%i, Monthly mean, # of av. rcrds=%i\n',imo,nrec);
  end

  fprintf('Saving MLD...\n');
  dmm = squeeze(MLD(mold).mld);
  MLD(mold).mld = dmm./nrec;

  fmat = sprintf('%smld_%i.mat',pthmat,yrold);
  fprintf('  ====  Saving %s\n',fmat);
  save(fmat,'MLD');
end    







