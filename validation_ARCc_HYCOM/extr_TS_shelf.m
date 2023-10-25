% Extract data for comparing with obs
% from HYCOM

addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
%addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_gom04/model_nest_ob;
startup

yr1=2008;
yr2=2008;
expt = 112;  % epxeriment with Greenland runoff and monthly Arctic rivers
regn = 'ARCc0.08';

%shlf='Beaufort';
shlf='Lapt';



rg=9806;  % convert pressure to depth, m
hg     = 2^100;  % "huge" in HYCOM used for land masking ONLY!
huge   = hg;     % 0-depth values are not = huge
hgg=1e20; 
btx = 'extr_TS_shelf.m';


pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthdata = '/Net/data2/ddmitry/Arctic_Data/BeaufortShelf/';
pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat2/',expt);

% -------------------------
% Get grid and bath:
% My bathymetry, Model bathymetry:
% -------------------------
ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

% Extract at observation locations: II, JJ
% Regions is contoured by IJ indices
% see beauf_shelf.m * Lapt_shelf.m
fobsmat=sprintf('%s%s_locations_%i-%i.mat',pthmat,shlf,yr1,yr2);
fprintf('Loading Obs. locations: %s\n',fobsmat);
OBS=load(fobsmat);
Iobs=OBS.II;
Jobs=OBS.JJ;
IndObs=sub2ind(size(HH),Jobs,Iobs);
IJ=OBS.IJ;

di=10;
[X,Y]=meshgrid([1:nn],[1:mm]);
IN=inpolygon(X,Y,IJ(:,1),IJ(:,2));
Ish0=find(IN==1 & HH<0);
Ish=Ish0(1:di:end);


f_plt=0;
if f_plt
figure(1); clf;
contour(HH,[0 0],'k');
hold on;
plot(IJ(:,1),IJ(:,2),'r-');
plot(X(Ish),Y(Ish),'g.');
plot(Iobs,Jobs,'m*');  % observations
end

id1=datenum(yr1,6,1)-datenum(yr1,1,1)+1;
id2=datenum(yr1,9,15)-datenum(yr1,1,1)+1;
dday=10;
for iyr=yr1:yr2
  yr=iyr;
  pthbin = sprintf('/nexsan/archive/%s_%3.3i/data/%i/',regn,expt,yr)
  if expt==112
    pthbin = sprintf('/nexsan/hycom/ARCc0.08_112/data/%i/',yr);
  end

  cc=0;
  SHLF=struct;
  SHLF.nm=shlf;
  for iday=id1:dday:id2
    dnmb=datenum(yr,1,1)+iday-1;
    fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
    finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);

    fprintf('Reading %s\n',datestr(dnmb));
    
    tic;
    [F,n,m,l] = read_hycom(fina,finb,'salin');
    F(F>hgg)=0;
    Ssh=squeeze(F(:,Ish)); % whole shelf
    SS=squeeze(F(:,IndObs));  % observation  locations

    [F,n,m,l] = read_hycom(fina,finb,'temp');
    F(F>hgg)=0;
    Tsh=squeeze(F(:,Ish));
    TT=squeeze(F(:,IndObs));

    [ZM0,ZZ0]=sub_zz_zm(fina,finb,HH,'f_btm',1);
    ZMsh=squeeze(ZM0(:,Ish));
%    ZZ=squeeze(ZZ0(:,Ish));
    ZM=squeeze(ZM0(:,IndObs));
%    ZZ=squeeze(ZZ0(:,IndObs));
%    [F,n,m,l] = read_hycom(fina,finb,'thknss');
%    F=F./rg;
%    F(F>hgg)=0;
%    dH=squeeze(F(1:Nlev,:,:)); 
    
    cc=cc+1;
    SHLF.TM(cc)=dnmb;
    SHLF.ZM(cc,:,:)=ZM;
    SHLF.Temp(cc,:,:)=TT;
    SHLF.Salin(cc,:,:)=SS;
    SHLF.Temp_sh(cc,:,:)=Tsh;
    SHLF.Salin_sh(cc,:,:)=Ssh;
    SHLF.ZM_sh(cc,:,:)=ZMsh;

    fprintf('Processed 1 reccord %4.1f min\n',toc/60);
  end
  SHLF.IJ_obs=IndObs;
  SHLF.IJ_shelf=Ish;
%  fmat=sprintf('%s%s-%3.3i_TSprf_shelf_%s_%i.mat',pthmat,regn,expt,shlf,yr)
  fmat=sprintf('%s%s-%3.3i_TSprf_obsloc_%s_%i.mat',pthmat,regn,expt,shlf,yr)
  fprintf('Saving %s\n',fmat);
  save(fmat,'SHLF');
end

