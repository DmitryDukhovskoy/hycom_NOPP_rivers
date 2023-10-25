% Preprare monthly river runoff data -> mat file
% with ready to use time series
% for creating HYCOM river fields 
%
% Using UCAR/NCAR CGD's  data set
% Dai and Trenberth data set
% http://www.cgd.ucar.edu/cas/catalog/surface/dai-runoff/index.html
% Many Arctic rivers have missing data 
% - interpolate in time
% ignore those with lots of missing data - these 
% are usually small rivers
% 
%
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom
addpath /usr/people/ddmitry/codes/MyMatlab/seawater
startup

clear 
close all

f_mat = 1; 

PTH.data    = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/force/rivers/';
PTH.dataout = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/force/riversNCAR/';
PTH.topo    = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
PTH.river   = '/Net/ocean/ddmitry/arctic_AOregimes/data/Greenland_rivers/';

pthncar = '/Net/ocean/ddmitry/UCAR_RIVERS/';
fnm = sprintf('%scoastal-stns-Vol-monthly.updated-Aug2014.nc',...
	      pthncar);
fmat = sprintf('%sncar_rivers_Arctic_1993-2016.mat',PTH.data);

Time = nc_varget(fnm,'time');
YRS  = floor(Time./100);
MM   = Time - YRS*100;
DD   = MM*0+1;
TM   = datenum(YRS,MM,DD);

% ARCc topo - check Topo version
% 11D - corrected version of T11, used in GOFS3.1
PTH.topo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
fltopo = sprintf('%sdepth_ARCc0.08_11.nc',PTH.topo);
HH  = nc_varget(fltopo,'Bathymetry');
LT = nc_varget(fltopo,'Latitude');
LN = nc_varget(fltopo,'Longitude');
[m,n]=size(HH);
[DX,DY] = sub_dx_dy(LN,LT);
ACell = DX.*DY;

% Find Arctic & Atlantic rivers & North Pacific
ocn = nc_varget(fnm,'ocn_name');
[a1,a2] = size(ocn);
rivname = nc_varget(fnm,'riv_name');
rlon = nc_varget(fnm,'lon_mou');
rlat = nc_varget(fnm,'lat_mou');
cc = 0;
rx = [];
for jk=1:a1
  ss = ocn(jk,1:3);
  srv = lower(rivname(jk,1:5));
%  if strncmp(srv,'neva',4); keyboard; end;
  if strncmp(ss,'ARC',3)
    cc = cc+1;
    rx(cc,1)=jk;
  elseif strncmp(ss,'ATL',3) & ...
	  rlat(jk)>50
    cc = cc+1;
    rx(cc,1)=jk;
  elseif strncmp(ss,'PAC',3) & ...
	  rlat(jk)>55
    cc = cc+1;
    rx(cc,1)=jk;
  end
      
end

% Khatanga, AWI/AOMIP climtology data
% for winter months - no observations
% just put some small values to avoid zeros
% Winter values are similar to Kolyma, which 
% has similar peak runoff
% Annual Nadym & Pur estimates from Zakharova's paper 2012
% Weights to convert annual to monthly rates
% derived from Taz (same location as Nadym and Pur)
% Anabar ~500 m3/s annual mean ~15.8 km3/yr
wght = [0.03, 0.02, 0.02, 0.02, 0.07, 0.37, 0.2, 0.07, ...
	0.07, 0.06, 0.04, 0.03];
Qkh=[100.0, 100.0, 100.0, 100.0, 100.0, 13378.5,...
     5766.38, 3241.54, 2825.38, 100.0, 100.0, 100.0]; % m3/s
Qpsn=[501.11, 501.11, 501.11, 501.11,501.11, 7516.72,...
      10022.3,2818.77, 3758.36, 4071.56, 1252.79, 689.03];
Qndm=wght*14.5*1e9/(3600*24*30.25); % km3/mo -> m3/s 
Qpur=wght*28*1e9/(3600*24*30.25); % km3/mo -> m3/s 
Qanb=wght*15.8*1e9/(3600*24*30.25); % km3/mo -> m3/s 

RVR.Nrivers  = length(rx);
RVR.Riv_name = rivname(rx,:);

it1  = find(TM == datenum(1993,1,1));
it2  = find(TM == datenum(2014,12,1));
Flow = nc_varget(fnm,'FLOW');
A    = Flow(it1:it2,rx);
[a1,a2] = size(A);
%keyboard
fprintf('Checking for completeness,IQ, and filling gaps...\n');
scs = 0;
for ik = 1:a2
  rnm = deblank(RVR.Riv_name(ik,:));
  Q = A(:,ik);
  Q(Q<=0)=nan;  % bad values

% Use AWI/AOMIP climatology for Khatanga  
  if strcmp(rnm,'Khatanga')
    Q(1:12)=Qkh;
    Q(13:24)=Qkh;
    Q(25:36)=Qkh;
  elseif strcmp(rnm,'Nadym')
    Q(1:12)=Qndm;
    Q(13:24)=Qndm;
    Q(25:36)=Qndm;
  elseif strcmp(rnm,'Pur')
    Q(1:12)=Qpur;
    Q(13:24)=Qpur;
    Q(25:36)=Qpur;
  elseif strcmp(rnm,'Anabar')
    Q(1:12)=Qanb;
    Q(13:24)=Qanb;
    Q(25:36)=Qanb;
% Add Pyasina River:
  elseif strcmp(rnm,'Poluy'); % instead of one of rivers to skip
    rnm = 'Pyasina            ';
    RVR.Riv_name(ik,1:length(rnm)) = rnm;
    rnm = deblank(rnm);
    Q = Q*nan;
    Q(1:12)=Qpsn;
    Q(13:24)=Qpsn;
    Q(25:36)=Qpsn;
    Q(37:48)=Qpsn;
%    keyboard
  end

%  if strcmp(rnm,'Severnaya Dvina');
%    keyboard;
%  end
  
% Skip tributaries and UK rivers
  rskp = sub_skip_rivers(rnm);
  if rskp, continue; end;

  Inan = find(isnan(Q));
  nmss = length(Inan);
  rmss = nmss/a1*100; 
  Rmss(ik,1) = rmss;
  fprintf('%i: %s, Missing data: %6.2f\n',...
	  ik,rnm,rmss);
  
  if rmss>94, % ~ 1year of data exists 
    fprintf('***  Skipping this river, too many missing data ...\n');
    continue; 
  end;
  
%  if rmss==0,
%    scs = scs+1;
%    Rname(scs,1:30) = RVR.Riv_name(ik,1:30);
%    Rindx(scs,1)    = rx(ik);
%    Qfxd(:,scs) = Q;
%    continue;
%  end
  
  if rmss>0 % fill gaps
    nyr = length(Q)/12;
    dmm = reshape(Q,[12,nyr]);
    Rmn = nanmean(dmm,2);

    jmnan = find(isnan(Rmn));
    if ~isempty(jmnan), 
      fprintf('***  Skipping this river, not enough data for monthly clim ...\n');
      continue; 
    end; % not enough data for monthly clim. 

    fprintf('%i:      Fixing gaps in data ...\n',ik);

    JJ = find(isnan(dmm));
    for jt=1:length(JJ)
      [j0,i0]=ind2sub(size(dmm),JJ(jt));
      dmm(j0,i0) = Rmn(j0);
    end
    qfxd = reshape(dmm,[12*nyr,1]);
    Q = qfxd;
  end
  
  if strcmp(rnm,'Lena'), % Split flow into 3 branches
    fprintf('Splitting Lena into 3 branches ...\n');
    scs=scs+1;
    Rname(scs,1:6) = 'Lena A';
    Rindx(scs,1)    = rx(ik);
    Qfxd(:,scs) = Q/3;
    
    scs=scs+1;
    Rname(scs,1:6) = 'Lena B';
    Rindx(scs,1)    = rx(ik);
    Qfxd(:,scs) = Q/3;
    
    scs=scs+1;
    Rname(scs,1:6) = 'Lena C';
    Rindx(scs,1)    = rx(ik);
    Qfxd(:,scs) = Q/3;
  elseif strcmp(rnm,'Mackenzie')
    fprintf('Splitting Mackenzie into 3 branches ...\n');
    scs=scs+1;
    Rname(scs,1:11) = 'Mackenzie A';
    Rindx(scs,1)    = rx(ik);
    Qfxd(:,scs) = Q/3;
    
    scs=scs+1;
    Rname(scs,1:11) = 'Mackenzie B';
    Rindx(scs,1)    = rx(ik);
    Qfxd(:,scs) = Q/3;
    
    scs=scs+1;
    Rname(scs,1:11) = 'Mackenzie C';
    Rindx(scs,1)    = rx(ik);
    Qfxd(:,scs) = Q/3;
  else
    scs = scs+1;
    Rname(scs,1:30) = RVR.Riv_name(ik,1:30);
    Rindx(scs,1)    = rx(ik);
    Qfxd(:,scs)=Q;
  end
  
end




fprintf('Total %i rivers identified & saved\n',scs);

%
% Extend Time series through Dec 2016:
% Using 2014 - last data in record, copy over
TMold = TM;
TM = TMold(it1:it2);
nrc = length(Qfxd);
Qfxd(nrc+1:nrc+12,:)=Qfxd(nrc-11:nrc,:);
nrc = length(Qfxd);
Qfxd(nrc+1:nrc+12,:)=Qfxd(nrc-11:nrc,:);

for yr=2015:2016
  for im=1:12
    TM(end+1)=datenum(yr,im,1);
  end
end

LAT = nc_varget(fnm,'lat_mou');
LON = nc_varget(fnm,'lon_mou');

clear RVR
RVR.Source = 'hycom_NOPP_rivers/prepare_rivers_arc08/fix_NCARrivers.m';
RVR.Info = 'River monthly flow, m3/sec, from NCAR Dai&Trenberth, extrapolated to2016';
RVR.TM = TM;
RVR.Riv_name = Rname;
RVR.Riv_indx = Rindx;
RVR.Qflow_m3_sec = Qfxd;
RVR.Lon = LON(Rindx);
RVR.Lat = LAT(Rindx);
Nrv = length(Rindx);

%RR = sub_find_river('Anabar',rivname,Flow,LAT,LON);
%RR = sub_find_river('Lena',rivname,Flow,LAT,LON);
%RR = sub_find_river('Pyasina',rivname,Flow,LAT,LON);
%RR = sub_find_river('Fish',rivname,Flow,LAT,LON);
%Qr = RR.Q;

% Adjust river mouth locations 
% to get correct positioning in ARCc
RVR = sub_correct_coordinates(RVR);

if f_mat>0
  fprintf('Saving rivers: %s\n',fmat);
  save(fmat,'RVR');
end

rnm = 'Stikine';
year = 2015;
mo = 7;
%[iC,jC,Qr] = sub_get_loc_rivername(rnm,LT,LN,HH,year,mo);
[TR,Qr] = sub_get_runoff_rname(rnm);
DR=datevec(TR);
Yr = [0:length(TR)-1]/12+DR(1);
figure(1); clf;
plot(Yr,Qr);
title(rnm);
set(gca,'xlim',[Yr(1) Yr(end)]);
btx = 'fix_NCARrivers.m';
bottom_text(btx,'pwd',1);

f_chck = 1;
if f_chck>0
  
  figure(10);
  clf;
  contour(HH,[0 0],'c');
  hold
  dmm=RVR.Lon;
  a1=length(dmm);
  for ik=1:a1
    x0=RVR.Lon(ik);
    y0=RVR.Lat(ik);
%    if y0<60, continue; end;
    D=distance_spheric_coord(LT,LN,y0,x0);
%    D=distance_spheric_coord(LAT,LON,y0,x0);
    [j0,i0]=find(D==min(min(D)));
    plot(i0,j0,'r*');
    nm = RVR.Riv_name(ik,1:8);
    text(i0,j0,nm);
  end
  btxt = 'fix_NCARrivers.m';
  bottom_text(btxt,'pwd',1);
  title('Actual locations of NCAR rivers, ARCc0.08');
  
end


  






