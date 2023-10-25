% Plot mean EKE
% from altimetry-derived U,V
% data provided by I. Bashmachnikov
%
% calculated in meanEKE.m
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

s_fig  = 0;

YR1=1993;
YR2=2016;

av=0; % =1 - average over specofoed years, =0 - plot individual years
plr =1;  % U from plr layer
rg = 9806;
rhow=1027; 
hgg=1e20;

regn = 'ARCc0.08';
expt = 110;
pthfig  = sprintf('/nexsan/people/ddmitry/Net_ocean/hycom/ARCc0.08/%3.3i/fig_meanUV/',expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat ='/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/data_mat/';


%fall=1; % overall mean EKE from mean <u'^2>,<v'2>, not same as mean[1993:2016]!

fprintf('\n Plotting altimetry mean EKE\n\n');

ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY] = sub_dx_dy(LON,LAT);

flalt1=sprintf('%sHHtransf_eke14_Greenl1993_2016S.mat',pthmat);
flalt2=sprintf('%sHHtransf_eke14_LsIrm1993_2016S.mat',pthmat);

fprintf('Loading %s\n',flalt1);
GR=load(flalt1);
fprintf('Loading %s\n',flalt2);
IR=load(flalt2);

A10=GR.log10eke;
A10(A10<-1e3)=NaN;
A=10.^(A10);
Alg=log(A);

B10=IR.log10eke;
B10(B10<-1e3)=NaN;
B=10.^(B10);
Blg=log(B);


% interpolate onto HYCOM grrid
fprintf('Interpolating to HYCOM region 1...\n');
[m1,n1]=size(Alg);
nmm=m1*n1;
kpp=0;
% subsample region to expedite calcu
sj1=380;
sj2=1100;
si1=800;
si2=1250;
slon=LON(sj1:sj2,si1:si2);
slat=LAT(sj1:sj2,si1:si2);

for jj=1:m1
  for ii=1:n1
    kpp=kpp+1;
    if mod(kpp,1000)==0,
      fprintf('===> %4.1f%%\n',kpp/nmm*100);
    end
    
    x0=GR.X(jj,ii);
    y0=GR.Y(jj,ii);
    dst=distance_spheric_coord(slat,slon,y0,x0);;
    [j1,i1]=find(dst==min(min(dst)));
    GR.I_hycom(jj,ii)=i1+si1-1;
    GR.J_hycom(jj,ii)=j1+sj1-1;
  end
end

% interpolate onto HYCOM grrid
fprintf('Interpolating to HYCOM region 2 ...\n');
[m2,n2]=size(Blg);
nmm=m2*n2;
kpp=0;
% subsample region to expedite calcu
bj1=1;
bj2=1000;
bi1=350;
bi2=1050;
slon=LON(bj1:bj2,bi1:bi2);
slat=LAT(bj1:bj2,bi1:bi2);

for jj=1:m2
  for ii=1:n2
    kpp=kpp+1;
    if mod(kpp,1000)==0,
      fprintf('===> %4.1f%%\n',kpp/nmm*100);
    end
    
    x0=IR.X(jj,ii);
    y0=IR.Y(jj,ii);
    dst=distance_spheric_coord(slat,slon,y0,x0);;
    [j1,i1]=find(dst==min(min(dst)));
    IR.I_hycom(jj,ii)=i1+bi1-1;
    IR.J_hycom(jj,ii)=j1+bj1-1;
  end
end


hps = [0.85 0.2 0.03 0.7];
c1=-7;
c2=-3;

% Greenland
xlim1 = si1;
xlim2 = si2;
ylim1 = sj1;
ylim2 = sj2;

stl = sprintf('Aaltimetry-derived, Annual Mean EKE, m2/s2, 1993-2016');
%  sub_plot_scalar_v2(S,nf,HH,xlim1,xlim2,ylim1,ylim2,LON,LAT,stl);
nf=1;
Lekem=Alg; % m2/s2
II=GR.I_hycom;
JJ=GR.J_hycom;
%pcolor(II,JJ,Alg); shading flat;
sub_plot_altmEKE(Lekem,II,JJ,nf,HH,...
		 xlim1,xlim2,ylim1,ylim2,...
		 LON,LAT,stl,...
		'c1',c1,'c2',c2,'clbpos',hps,'cmp',5);
%contour(ekem,[10:10:100],'k');

btx='plot_EKEaltimetry.m';
bottom_text(btx,'pwd',1,'position',[0.05 0.08 0.4 0.04]);

if s_fig==1
  fgnm=sprintf('%saltimetry_meanEKEreg1_1993-2016',pthfig);
  fprintf('Saving %s\n',fgnm);
  print('-dpng','-r200',fgnm);
end  


% Irminger
xlim1 = bi1;
xlim2 = bi2;
ylim1 = bj1;
ylim2 = bj2;

stl = sprintf('Aaltimetry, Annual Mean EKE, m2/s2, 1993-2016');
%  sub_plot_scalar_v2(S,nf,HH,xlim1,xlim2,ylim1,ylim2,LON,LAT,stl);

nf=2;
Lekem=Blg; % m2/s2
II=IR.I_hycom;
JJ=IR.J_hycom;
%pcolor(II,JJ,Alg); shading flat;
%clbtck=[c1:0.6:c2];
%u=exp(clbtck)*10000; %cm2/s2
%for kk=1:length(clbtck)
%  clblbl{kk}=sprintf('%5.1f',u(kk));
%end

sub_plot_altmEKE(Lekem,II,JJ,nf,HH,...
		 xlim1,xlim2,ylim1,ylim2,...
		 LON,LAT,stl,...
		'c1',c1,'c2',c2,'clbpos',hps,'cmp',5);
%contour(ekem,[10:10:100],'k');

btx='plot_EKEaltimetry.m';
bottom_text(btx,'pwd',1,'position',[0.05 0.08 0.4 0.04]);

if s_fig==1
  fgnm=sprintf('%saltimetry_meanEKEreg2_1993-2016',pthfig);
  fprintf('Saving %s\n',fgnm);
  print('-dpng','-r200',fgnm);
end  



