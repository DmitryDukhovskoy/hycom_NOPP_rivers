% Plot mean SSH from GLobal HYCOM reanalysis
% GLBb data were remapped onto ARCc grid
% These files include GLBb reanalysis (for 1993-2012)
% and analysis (expt 90.[01])
% see: ~ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers
%  /REMAP_ARCc/remap_gridGLB2ARC_archv/remap_annualNyrs_analysisGLBa.csh
% integrating to the depth of Sref
% /nexsan/GLBa0.08/expt_90.9/data
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

s_fig  = 1;


rg=9806;  % convert pressure to depth, m
Sref=35; % N.Atl. is too saline
TV = '07';

pthARC = '/Net/mars/ddmitry/hycom/GLBb2ARCc0.08/mean/';
pthmat = '/Net/mars/ddmitry/hycom/GLBb0.08/data_mat/';
pthfig = '/Net/mars/ddmitry/hycom/GLBb0.08/fig_ssh/';
%ftopo = sprintf('%sGLBb_T07_subset_Natl.mat',pthmat);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
%pthglb = '/nexsan/GLBb0.08/GLBb0.08_191/topo/';

YRPLT=[];
cnc=0;
for yr=1993:2015
  cnc=cnc+1;
  YRPLT(cnc,1)=yr;
end

ftopo = sprintf('%sdepth_ARCc0.08_%s.nc',pthtopo,TV); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');

IND = sub_get_ARCgrid;
ind1=IND.i1;
ind2=IND.i2;
jnd1=IND.j1;
jnd2=IND.j2;

HH  = HH(jnd1:jnd2,ind1:ind2);
LON = LON(jnd1:jnd2,ind1:ind2);
LAT = LAT(jnd1:jnd2,ind1:ind2);
[mm,nn]=size(HH);


% Mask of region of interest:
LMSK=HH;
LMSK(LMSK>=-5)=0;
LMSK(LMSK<0)=1;
IDEEP=find(LMSK>0);
INAN =find(LMSK==0);

hmin = 0;
ARC = sub_arctic_domain_glb008(HH,hmin);

%Hmsk = HH*0;
%Hmsk(HH>0)=1;
%cmm = [1,1,1; 0.6,0.6,0.6];


np=length(YRPLT);
iy1=1;
cc=iy1-1;
for iyr=iy1:np
  tic;
  yr=YRPLT(iyr);
  fprintf('Reading SSH %i\n',yr);
  
  pthbin=pthARC;
  if yr<1995
    expt=190;
  elseif yr>=1995 & yr <2013
    expt=191;
  elseif yr==2013
    expt=910; % 
  elseif yr==2014 | yr==2015
    expt=911; % 
  end
  
  fnm=sprintf('%s%3.3i_archMN_GLBb2ARCc.%4.4i_01_%4.4i_12',pthbin,expt,yr,yr);
  fina=[fnm,'.a'];
  finb=[fnm,'.b'];
  
  
  fld='srfhgt';
  [F,nn,mm,ll] = read_hycom(fina,finb,fld);
  toc;
  F(F>1e6)=nan;
  F=squeeze(F/9.806);
  E=F(jnd1:jnd2,ind1:ind2); % subset

  Inan = ARC.Inan;
  E(Inan)=nan;
  E=E-nanmean(nanmean(E));

  nf = -1;
  sub_plot_ssh(E,HH,LON,LAT,nf);

  ttl=sprintf('Annual SSH, Global HYCOM+NCODA, %4.4i ',yr);
  title(ttl);
  btx = 'plot_meanSSH.m';
  bottom_text(btx,'pwd',1);
  
  if s_fig==1
    fgnm=sprintf('%sglb008_meanSSH%4.4i',pthfig,yr);
    fprintf('Saving %s\n',fgnm);
    print('-dpng','-r250',fgnm);
  end
  

end

