% Interpolate HDEM fields onto HYCOM grid 0.08
% Nearest point interpolation - fastest for plotting only
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08/;
addpath /usr/people/ddmitry/codes/MyMatlab/;
startup;

clear
close all

f_smat = 1;

pthmat  = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/data_mat/';
pthin   = '/Net/data/GDEM4/';  % climatology data with no land
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
ffout   = sprintf('%sGDEM_gridindx_arc008.mat',pthmat);

ftopo = sprintf('%sdepth_ARCc0.08_11.nc',pthtopo);
HT  = nc_varget(ftopo,'Bathymetry');
XT = nc_varget(ftopo,'Longitude');
YT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(HT);
[IT,JT] = meshgrid((1:nn),(1:mm));

%[DX,DY]=sub_dx_dy(XT,YT);

fprintf('loading HYCOM indces for GDEM %s\n',ffout);
load(ffout);
IWP = INDX.IWP;
JWP = INDX.JWP;


mo=8;

tnc = sprintf('%sptgdemv4f%2.2i.nc4',pthin,mo);
snc = sprintf('%ssgdemv4f%2.2i.nc4',pthin,mo);
zz = nc_varget(tnc,'Depth');
nz = length(zz);
T = nc_varget(tnc,'Potential_Temperature');
S = nc_varget(snc,'salinity');

fprintf('Start interploation GDEM--> 0.08 HYCOM Arctic, month=%i\n',mo);

Ib = find(IWP>0 & JWP>0);
Igdem = round(IWP(Ib));
Jgdem = round(JWP(Ib));
TT = squeeze(T(:,Ib));
SS = squeeze(S(:,Ib));

Lmsk = HT*0;
Lmsk(HT<0) = 1;
Iocn = find(Lmsk==1);

%Sav = HT*nan;
%Tav = HT*nan;
for ii=1:length(Iocn)
  if mod(ii,100000)==0,
    fprintf('Interpolating %5.1f%% ...\n',ii/length(Iocn)*100);
  end

  lix = Iocn(ii);
  [j0,i0] = ind2sub(size(HT),lix);

  d0=sqrt((Igdem-i0).^2+(Jgdem-j0).^2);
  imin = find(d0==min(d0),1);
  if d0(imin)>10;
    Sav(1:nz,j0,i0) = nan;
    continue;
  end
 
  Sav(:,j0,i0) = SS(:,imin);
  Tav(:,j0,i0) = TT(:,imin);

 
end 

fmatT=sprintf('%sTave_GDEM_%2.2i.mat',pthmat,mo)
fmatS=sprintf('%sSave_GDEM_%2.2i.mat',pthmat,mo)

fprintf('Saving %s\n',fmatT);
save(fmatT,'-v7.3','Tav','zz');

fprintf('Saving %s\n',fmatS);
save(fmatS,'-v7.3','Sav','zz');


