% Plot monthly mean UV fields
% derived in mnthly_arc08_UV50m.m
% The velocity is averaged over layers zz1-zz2
% 
% <EKE> = 1/2*(u^2+v^2)-1/2(<u^2>+<v^2>)=
% ...=1/2*(<u'^2>+<v'^2>)  (1)
% where <*> is averaging over long time scale
%
% if instanteneous output are not available
% <EKE>=1/2[(u^2-<u^2>)+(v^2-<v^2>)]  (2)
% 
% I am calculating directly EKE using (1)
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

s_fig  = 0;
s_mat = 1; % save overall mean; =0 - load saved and plot

dday=5;
rg = 9806;
rhow=1027; 
hgg=1e20;

zz1=0;
%zz2=-50;
zz2 = -15;


regn = 'ARCc0.08';
expt = 110;
pthfig  = '/nexsan/people/ddmitry/Net_ocean/hycom/ARCc0.08/110/fig_meanUV/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat =sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
pthmat2 =sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat2/',expt);

YR1=2009;
YR2=YR1;
YRPLT = [YR1:YR2];
np = length(YRPLT);

fprintf('Calculating mean EKE, %i-%i, %i - %i m\n',YR1,YR2,abs(zz1),abs(zz2));
fprintf('EKE = 1/2*<(u'')^2+(v'')^2> \n');


ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY] = sub_dx_dy(LON,LAT);

%U2MN=HH*0; % all-years mean
%V2MN=HH*0;
%ccT=0;
for ikk=1:np
  iyr = YRPLT(ikk);
%  fmat = sprintf('%smnthUV_lr%2.2i_%i.mat',pthmat,plr,iyr);
  fmat = sprintf('%smnthUV_%4.4i-%4.4i_%i.mat',...
		 pthmat2,abs(zz1),abs(zz2),iyr);
  fprintf('Loading %s\n',fmat);
  load(fmat);

  fout = sprintf('%sarc08_%i_annualEKE_%4.4i-%4.4i_%i.mat',...
	       pthmat2,expt,abs(zz1),abs(zz2),iyr);
  
% Annual means:
  cc = 0;
  u2mn=zeros(mm,nn);  % 1-yr mean
  v2mn=zeros(mm,nn);
  for iday=1:dday:365
    tic;
    dnmb = datenum(iyr,1,1)+iday-1;
    DV   = datevec(dnmb);
    imo = DV(2);

    pthbin = sprintf('/nexsan/archive/ARCc0.08_%3.3i/data/%4.4i/',expt,iyr);
    fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,iyr,iday);
    finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,iyr,iday);

    if ~exist(fina,'file') | ~exist(finb,'file')
      fprintf('Not found %s or %s\n\n',fina,finb);
      continue;
    end

% Mean fields:    
    Umn = meanUV(imo).U;
    Vmn = meanUV(imo).V;

    fprintf('Reading %4.4i/%2.2i/%2.2i: %s\n',DV(1:3),fina);

    [F,n,m,nlr] = read_hycom(fina,finb,'u-vel.');
    F=squeeze(F);
    F(F>hgg)=nan;
    U=F;
%    F=squeeze(F(1:plr,:,:));
%    U=squeeze(nanmean(F,1));

%    [F,n,m,nlr] = read_hycom(fina,finb,'v-vel.','r_layer',plr);
    [F,n,m,nlr] = read_hycom(fina,finb,'v-vel.');
    F=squeeze(F);
    F(F>hgg)=nan;
    V=F;
%    F=squeeze(F(1:plr,:,:));
%    V=squeeze(nanmean(F,1));

% Average over specified depths:    
    fprintf('Getting layer depths ...\n');
    [ZM,ZZ] = sub_zz_zm(fina, finb,HH,'f_btm',1);

    dz=abs(zz2-zz1);
    zbtm=-9000;
    Uav = sub_zLayer_average(HH,ZZ,U,zbtm,zz1,zz2); % depth-average field 
    Vav = sub_zLayer_average(HH,ZZ,V,zbtm,zz1,zz2); % depth-average field 

% Need to collocate U, V and dH
%  [ll,mm,nn]=size(U);
    Uclc=Uav*0;
    Vclc=Vav*0;
    for ii=2:nn-1
      u1 = Uav(:,ii);
      u2 = Uav(:,ii+1);
      uu = 0.5*(u1+u2);
      Uclc(:,ii) = uu; % collocated U in the layer
    end
    for jj=2:mm-1
      v1 = Vav(jj,:);
      v2 = Vav(jj+1,:);
      vv = 0.5*(v1+v2);
      Vclc(jj,:)=vv;
    end
    dU = Uclc-Umn;
    dV = Vclc-Vmn;

    Uprt2 = dU.^2;
    Vprt2 = dV.^2;
    
% !!!  eke = 1/2*rhow*(dU.^2+dV.^2); This is wrong formula !!!
% should be 1/2*rhow*(mean(dU^2)+mean(dV^2)) !!!
    cc=cc+1;
    u2mn=u2mn+Uprt2;
    v2mn=v2mn+Vprt2;

    % Diagnostics
    umm=u2mn/cc; % <u'^2>
    vmm=v2mn/cc; % <v'^2>
    eke=1/2*(umm+vmm);
    fprintf('Mean EKE=%8.4f m2/s2, Max EKE=%8.4f m2/s2\n',...
	    nanmean(nanmean(eke)),max(max(eke)));
    fprintf('++++++>  Processed 1 record %8.5f min\n\n',toc/60);

  end  % days

  u2mn=u2mn/cc; % <u'^2>
  v2mn=v2mn/cc; % <v'^2>
  ekem=1/2*(u2mn+v2mn);
  
  EKE.Title=sprintf('Annual mean EKE from HYCOM-CICE 0.08-110, %i-%i m averg U',...
             abs(zz1),abs(zz2));
  EKE.Formula='<eke>=1/2*(<u''^2>+<v''^2>)';
  EKE.Year=iyr;
  EKE.Number_records_averaged=cc;
  EKE.EKE_m2_s2 = ekem;  % annual mean

  fprintf('Saving mean %s\n',fout);
  save(fout,'EKE');

end

%U2MN=U2MN/ccT;
%V2MN=V2MN/ccT;
%ekeM=1/2*(U2MN+V2MN);
%EKE(ikk+1).Title='OVERALL 1993-2016 mean EKE from HYCOM-CICE 0.08-110';
%EKE(ikk+1).Formula='<eke>=1/2*(<u''^2>+<v''^2>)';
%EKE(ikk+1).Year=iyr;
%EKE(ikk+1).Number_records_averaged=ccT;
%EKE(ikk+1).EKE_m2_s2 = ekeM;  % total mean
