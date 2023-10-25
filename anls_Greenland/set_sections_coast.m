% For Greenland tracer analysis
% and estimating across-flow eddy fluxes 
% define segments
% along Gr. coast

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

regn = 'ARCc0.08';
expt = '110';
%TV   = 11; % topo
TV   = 11; % topo
pthfig  = sprintf('/Net/mars/ddmitry/hycom/ARCc0.08/%s/fig_domain/',expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = '/Net/mars/ddmitry/hycom/ARCc0.08/data_mat/';
%pthmat  = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/tmp_data/';

s_fig=0;
s_mat=1;

txtb='ddmitry@mars: hycom_NOPP_rivers/anls_Greenland/set_sections_coast.m';

ftopo = sprintf('%s/depth_%s_%i.nc',pthtopo,regn,TV); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

domname = 'NorthAtl';
%IND = smaller_domain_indices('Green');
IND = smaller_domain_indices(domname);
%IND = smaller_domain_indices('subpolgyre');
%IND = [];


if isempty(IND) % whole domain
  IND.i1 = 1;
  IND.i2 = nn;
  IND.j1 = 1;
  IND.j2 = mm;
  IND.dj = mm;
  IND.di = nn;
end
  
inc1=IND.i1;
inc2=IND.i2;
jnc1=IND.j1;
jnc2=IND.j2;
djnc=IND.dj;
dinc=IND.di;


HH=HH(jnc1:jnc2,inc1:inc2);
LON=LON(jnc1:jnc2,inc1:inc2);
LAT=LAT(jnc1:jnc2,inc1:inc2);
[mm,nn]=size(LON);


fn=1;
sub_plot_bath(HH,LON,LAT,fn,domname);

% For defining norms and sign convention
% need direction to the coast
IJC=[   402   955
   308   907
   290   861
   303   795
   272   686
   282   620
   249   560
   215   526
   215   383
   261   362
   306   439
   321   475
   462   550
   532   839
   475   917
   422   930];
IJC(end+1,:)=IJC(1,:);
nsgm = length(IJC);
ii=0;
for isgm=1:nsgm-1
  is=IJC(isgm,1);
  ie=IJC(isgm+1,1);
  js=IJC(isgm,2);
  je=IJC(isgm+1,2);
  [II,JJ]=sub_xsct_indx(is,js,ie,je);
  plot(II,JJ,'c.');
  I=sub2ind(size(HH),JJ,II);

  ii=ii+1;
  CST(ii).Name=sprintf('CoastLine_%2.2i',ii);
  CST(ii).I=II;
  CST(ii).J=JJ;
  CST(ii).LON=LON(I);
  CST(ii).LAT=LAT(I);
end;
% Combine into 1 contour:
IJC=[];
for ii=1:nsgm-1
  I=CST(ii).I';
  J=CST(ii).J';
  dmm=[I,J];
  IJC=[IJC;dmm];
end


% Near the coast:
IJ=[   281   798
   237   678
   190   552
   188   436
   203   359
   258   312
   321   423
   457   500
   529   545
   535   629
   571   681
   578   834];

ii=0;
nsgm = length(IJ);
for isgm=1:nsgm-1
  is=IJ(isgm,1);
  ie=IJ(isgm+1,1);
  js=IJ(isgm,2);
  je=IJ(isgm+1,2);
  [II,JJ]=sub_xsct_indx(is,js,ie,je);
  plot(II,JJ,'b.-');
  I=sub2ind(size(HH),JJ,II);

  [nx,ny]=sub_find_norm(II,JJ,IJC); 
  nx=sign(nx);
  ny=sign(ny);
  
  ii=ii+1;
  SEGM(ii).Name=sprintf('Contour1_%2.2i',ii);
  SEGM(ii).I=II;
  SEGM(ii).J=JJ;
  SEGM(ii).LON=LON(I);
  SEGM(ii).LAT=LAT(I);
  SEGM(ii).Norm=[nx,ny]; % local normal, determines + flux [X,Y]
end;


% Off the coast:
IJ=[   261   808
   213   689
   174   558
   172   436
   188   350
   258   294
   339   409
   467   480
   550   541
   564   621
   604   677
   613   837];
%plot(IJ(:,1),IJ(:,2),'r.-');

nsgm = length(IJ);
for isgm=1:nsgm-1
  is=IJ(isgm,1);
  ie=IJ(isgm+1,1);
  js=IJ(isgm,2);
  je=IJ(isgm+1,2);
  [II,JJ]=sub_xsct_indx(is,js,ie,je);
  plot(II,JJ,'r.-');
  I=sub2ind(size(HH),JJ,II);

  [nx,ny]=sub_find_norm(II,JJ,IJC); 
  nx=sign(nx);
  ny=sign(ny);
%  if isgm==8, keyboard; end;
  ii=ii+1;
  SEGM(ii).Name=sprintf('Contour2_%2.2i',ii); % 2nd contour
  SEGM(ii).I=II;
  SEGM(ii).J=JJ;
  SEGM(ii).LON=LON(I);
  SEGM(ii).LAT=LAT(I);
  SEGM(ii).Norm=[nx,ny]; % local normal, determines + flux [X,Y]
end;




set(gca,'xlim',[50 800],'ylim',[250 1000]);

bottom_text(txtb);
if s_fig==1
  fgnm=sprintf('%ssections_GrCoast',pthfig);
  fprintf('Saving fig %s\n',fgnm);
  print('-dpng','-r200',fgnm);
end


if s_mat==1
  fmat=sprintf('%sGrCoast_sections',pthmat);
  fprintf('Saving %s\n',fmat);
  save(fmat,'SEGM','IND','IJC');
end

