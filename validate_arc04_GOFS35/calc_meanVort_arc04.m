% Calculate area-mean vorticity
% ocean currents
% using circulation theorem
%
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear


YR1=2020;
YR2=YR1;
dday=5;

RR = 200; % ring size (diamter), km 
Ilr = 5;  % Z layer where vort is calcualted

%regn = 'natl'; % Natl region
%regn = 'arctA'; % 
%regn = 'arctB'; % 
regn = 'ArctOc'; % 


% This code is for archv output fields
% CHeck if archm - do not add U/V barotropic !!! - use extr_*mean_ 
s_mat=1;  % =2 - load saved and start from the last record

ixx    = 9; % experiment name and dir - check with EXPT - expt 023
%ixx    = 6;  % expt 022 original 
EXPT   = sub_cice_experiments;
expt   = EXPT(ixx).Nmb;
texpt  = EXPT(ixx).cice_opt; % CICE options for sens. experiments
res    = EXPT(ixx).res;


%Zmn = 100; % average over the top Zmn m
rg  = 9806;
hgg = 1e20; % 


pthmat = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/%3.3i/strait_fluxes/',expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';

YRPLT=[];
cc=0;

%if expt=='061', yend=2008; end;
dd = 5;
for iyr=YR1:YR2
  nday=datenum(iyr,12,31)-datenum(iyr,1,1)+1;
  for iday=2:dd:nday
    cc=cc+1;
    YRPLT(cc,1)=iyr;
    YRPLT(cc,2)=iday;
    j1d  = datenum(iyr,1,1);
    dnmb = j1d+iday-1;
    YRPLT(cc,3)=dnmb;
  end
end

% Start from:
ftopo = sprintf('%sdepth_ARCc0.04_17DD.nc',pthtopo); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');

[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);
[II,JJ] = meshgrid([1:nn],[1:mm]);


switch(regn)
 case('ArctOc')
  IJr = [ 720        3734
         760        3167
         882        2933
        1060        2750
        1161        2541
        1273        2425
        1415        2333
        1491        2274
        1649        2124
        1806        2041
        1892        1957
        2167        1949
        2314        2233
        2405        2575
        2446        3067
        2217        3751
        1390        3868];

end

xl1 = min(IJr(:,1));
xl2 = max(IJr(:,1));
yl1 = min(IJr(:,2));
yl2 = max(IJr(:,2));

% Subset grid points:
dii = 6;
HHp=HH*0;
HHp(yl1:dii:yl2,xl1:dii:xl2)=HH(yl1:dii:yl2,xl1:dii:xl2);
[IIp,JJp] = meshgrid([xl1:dii:xl2],[yl1:dii:yl2]);

hmin = -500;

% Define points in the region
inp = inpolygon(II,JJ,IJr(:,1),IJr(:,2));
%inp_sub = inpolygon(IIp,JJp,IJr(:,1),IJr(:,2));

% Define deep ocean points only
Ioc = find(inp==1 & HHp<hmin);
%Inc = find(inp~=1 & HH>=hmin);

CRL = [];
CINDX = [];
cc = 0;
for iyr = YR1:YR2
  yr = iyr;

  dE=datenum(yr,12,31);
  dJ1=datenum(yr,1,1);
  ndays=dE-dJ1+1;

	fmat = sprintf('%s%3.3i_meanVort_%s_%4.4i.mat',...
					 pthmat,expt,regn,iyr);
	fprintf('Data will be saved -> %s\n',fmat);

  dlast=0;
  if s_mat == 2
    fprintf('Loading saved %s start from last recrod\n',fmatout);
    load(fmat);
    dlast = SCT(1).Time(end);
    fprintf('Last saved record %s\n\n',datestr(dlast));
  end


  iday1=1;
  iday2=ndays;

  for iday=iday1:dday:iday2
    iday0=iday;
    dnmb=datenum(dJ1)+iday0-1;
    DV   = datevec(dnmb);

%    pthbin=sprintf('/nexsan/hycom/ARCc0.04_%3.3i/data/%4.4i/',expt,yr);
%    pthbin=sprintf('/nexsan/people/ddmitry/hycom/ARCc0.04_022/data/%4.4i/',yr)
    pthbin=sprintf('/nexsan/people/ddmitry/hycom/ARCc0.04_%3.3i/data/%4.4i_%s/',...
                    expt,yr,texpt);

    fina = sprintf('%s%3.3i_archv.%4.4i_%3.3i_00.a',pthbin,expt,yr,iday);
    finb = sprintf('%s%3.3i_archv.%4.4i_%3.3i_00.b',pthbin,expt,yr,iday);

    if yr==2017 & DV(2)<4
      pthbin='/nexsan/people/ddmitry/hycom/ARCc0.04_022/data/2017_BL99Tfrz/';
      fina = sprintf('%s022_archv.%4.4i_%3.3i_00.a',pthbin,yr,iday);
      finb = sprintf('%s022_archv.%4.4i_%3.3i_00.b',pthbin,yr,iday);
    end


%
% Search for close days:
		for ik=1:dday
			iday = iday-1;
			dnmb = datenum(dJ1)+iday-1;
			DV   = datevec(dnmb);
			fina = sprintf('%s%3.3i_archv.%4.4i_%3.3i_00.a',pthbin,expt,yr,iday);
			finb = sprintf('%s%3.3i_archv.%4.4i_%3.3i_00.b',pthbin,expt,yr,iday);
			if exist(fina,'file');
				fprintf(' Found closest file: %s\n',fina);
				break;
			end;
		end
% If went back to iday0 - skip
		if ik==dday & ~exist(fina,'file');
			continue;
		end

    cc=cc+1;
    if dnmb<=dlast
      fprintf('Record exist %s, skipping ...\n',datestr(dnmb));
      continue;
    end
    fprintf('Reading %4.4i/%2.2i/%2.2i: %s\n',DV(1:3),fina);

    tic;
    [F,n,m,nlr] = read_hycom(fina,finb,'u-vel.','r_layer',Ilr);
    F(F>hgg)=0;
    U=squeeze(F);

    [F,n,m,llr] = read_hycom(fina,finb,'u_btrop');
    F(F>hgg)=0;
    Ub=squeeze(F);

    U=U+Ub;

% These are NOT total velocities    - archv
% barotropic needed
    [F,n,m,nlr] = read_hycom(fina,finb,'v-vel.','r_layer',Ilr);
    F(F>hgg)=0;
    V=squeeze(F);

    [F,n,m,llr] = read_hycom(fina,finb,'v_btrop');
    F(F>hgg)=0;
    Vb=squeeze(F);

    V=V+Vb;

% Mean curl over Nordic Seas getpts
    [VORT,CINDX] = vort_areamean_hycom(LON,LAT,U,V,Ioc,II,JJ,RR,CINDX);
    VRT = VORT.VRT;

    dmm = VORT.VRT;
    MVRT.VRT(cc,:) = dmm;
    MVRT.Units = 's-1';
    MVRT.TM(cc,1) = dnmb;
    MVRT.Ioc = Ioc; 

    fprintf('Processing 1 time rec: %8.3f sec\n',toc);

    if mod(cc,10)==0
      fprintf('  ====  Saving %s\n',fmat);
      save(fmat,'MVRT');
    end

  end
end

fprintf('  =====  End  ====  Saving %s\n',fmat);
save(fmat,'MVRT');

% Check
f_chck=0;
if f_chck==1
	Ioc = MVRT.Ioc;
	ik = 1;
	AA = HH*0;
	AA(HH>=0)=nan;
	AA(Ioc) = MVRT.VRT(ik,:); 
	clear Vrt
	[ms,ns]=size(IIp);
	for isb=1:ns
		for jsb=1:ms
			isb0=IIp(jsb,isb);
			jsb0=JJp(jsb,isb);
			Vrt(jsb,isb)=AA(jsb0,isb0);
		end
	end

end




