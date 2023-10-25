% annual mean S in North Atlantic
% Updated verions with bias correction
% is used
% When using the data set in a paper, please quote 
% the version number (EN.4.2.1), the date the data were downloaded 
% and cite the following:
%
% http://hadobs.metoffice.com/en4/index.html
%The EN4 dataset consists of two products:
%Observed subsurface ocean temperature and 
% salinity profiles with data quality information, and,
%Objective analyses formed from the profile data with uncertainty estimates.
%Data are available from 1900 to the present and 
% there are separate files for each month.
% Please read 'Good, S. A., M. J. Martin and N. A. Rayner, 2013. EN4: 
% quality controlled ocean temperature and salinity profiles and 
% monthly objective analyses with uncertainty estimates, 
% Journal of Geophysical Research: Oceans, 118, 6704-6716, 
% doi:10.1002/2013JC009067' for details of how the dataset was constructed.
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

s_fig = 0; 
s_mat = 1; % =0 - load existing mat file
           % =1 calculate annual S


%pthdat = '/Net/yucatan/tachanat/ocean_analysis/EN4/EN4_extract/';
pthdat = '/Net/kronos/ddmitry/EN4/';
%pthmat = '/Net/ocean/ddmitry/vector_winds/dataM/'; 
pthmat = '/nexsan/people/ddmitry/data_mat/';
pthmat2= '/Net/mars/ddmitry/hycom/GLBb0.08/data_mat/';
pthfig = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/fig_EN4/';
en4v   = 'EN.4.2.0.f.analysis.g10'; % EN4 version
txtb   = 'annualS_subarctic.m';

LR=[0, -51;...
    -51,-151;...
    -151,-301];
nlrs=size(LR,1);

ftopo = sprintf('%sGLBb_T07_subset_Natl.mat',pthmat2);
fprintf('Loading topo %s\n',ftopo);
TH = load(ftopo);
HH = TH.HH;
LONH = TH.LON;
LATH = TH.LAT;

%Sref = 35.5;
%Sref = 35; 
yr1=1990;
yr2=2016;
fmat = sprintf('%sannualS_subarctic_EN4.mat',pthmat);


fnm = sprintf('%s%s.201302.nc',pthdat,en4v);

S = double(nc_varget(fnm,'salinity'));
S = squeeze(S(1,1,:,:));

% Find sections:
LON = nc_varget(fnm,'lon');
LAT = nc_varget(fnm,'lat');
ZM  = -1*nc_varget(fnm,'depth');
lz = length(ZM);
dzm = abs(diff(ZM));
dzm(lz) = dzm(lz-1);

clear ZZ
ZZ(1,1) = 0;
for kk=1:lz
  ZZ(kk+1) = -(abs(ZM(kk))+abs(0.5*dzm(kk)));
end
ZZ=ZZ(:);
DZ=abs(diff(ZZ)); % layer thicknesses

% Find layers;
for ilv=1:nlrs
  z1=LR(ilv,1);
  z2=LR(ilv,2);
  ik1=max(find(ZZ>=z1));
  ik2=max(find(ZZ>=z2));
  iLR(ilv,1)=ik1;
  iLR(ilv,2)=ik2;
end

I = find(LON>180);
LON(I)=LON(I)-360;

%FWC = sub_regions;
%nR = length(FWC);

%keyboard


[LN,LT] = meshgrid(LON,LAT);
[DX,DY] = sub_dx_dy(LN,LT);
Acell = DX.*DY*1e-6; % km2


cc=0;
YRPLT=[];
for ii=yr1:yr2
    cc=cc+1;
    YRPLT(cc,1)=ii;
end
nyrs=cc;


for ik=1:nyrs
  YR=YRPLT(ik,1);

  cc=0;
  for IM=1:12
    dnmb = datenum(YR,IM,1);
    fnm = sprintf('%s%s.%4.4i%2.2i.nc',pthdat,en4v,YR,IM);

    fprintf('EN4: Reading %i/%i, %s\n',YR,IM,fnm);

    SAL = double(nc_varget(fnm,'salinity'));
%    TT = double(nc_varget(fnm,'temperature'))-273.15; % K -> C
    cc=cc+1;

% North of 30N    
    SAL = squeeze(SAL(1,:,120:end,:));
%    TT = squeeze(TT(1,:,:,:));
    [ll,mm,nn]=size(SAL);

    if IM==1,
      sumS=zeros(nlrs,mm,nn);
    end
    
    for ilv=1:nlrs
      k1=iLR(ilv,1);
      k2=iLR(ilv,2);
      
      smm = zeros(mm,nn);
      sdz=0;
      for kk=k1:k2-1
	S=squeeze(SAL(kk,:,:));
	dz=abs(DZ(kk));
	sdz=sdz+dz;
	smm=smm+dz*S;
      end
      
      sumS(ilv,:,:)=squeeze(sumS(ilv,:,:))+smm./sdz; % sum by month, depth-average S
    end
%keyboard

  end;  %months
  
  sumS=sumS/cc;
  
  fprintf(' Saving annual S, # mo=%i, maxS=%6.2f minS=%6.2f\n',...
	  cc,max(max(sumS(1,:,:))), min(min(sumS(1,:,:))));
  
  SANN(ik).vert_layers = LR;
  SANN(ik).year        = YR;
  SANN(ik).annualS     = sumS;
  SANN(ik).LON         = LON;
  SANN(ik).LAT         = LAT(120:end);
  
  if s_mat==1
    fprintf('Saving %s\n',fmat);
    save(fmat,'SANN');
  end

end  % years

  

  
