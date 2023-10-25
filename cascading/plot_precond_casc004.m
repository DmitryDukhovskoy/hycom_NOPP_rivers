% Plot locations and probability
% of surf-2-bottom
% convection
% see: precondition_casc.m - saved by winter/year (Oct-Mar yr1/yr2)=yr1
%
% To be able to be cascaded downslope in the benthic layer, 
% the mixed layer depth (MLD) on the shelf-slope must reach
% the bottom
%
% expt_110 - no Greenland runoff
% expt 112 - with Greenland runoff
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

% if yr2>yr1 - combine statistics over specified years
yr1=2005;
yr2=2005;

regn = 'ARCc0.04';
%expt = 011; % no Greenland runoff  
%expt = 012;  % Greenland runoff

rg = 9806;
% Thershold dRho
dRho0=1e-4*1023;

fprintf('Plotting preconditioning casc: %s-%3.3i %i-%i\n\n',...
	regn,expt,yr1,yr2);

pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/topo_grid/';
pthmat=sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/%s/%3.3i/data_mat/',regn,expt);
btx = 'plot_precond_casc004.m';



ftopo = sprintf('%s/depth_%s_17DD.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
%[XX,YY]=meshgrid((1:nn),(1:mm));
cnc=0;
CONV=HH*0;
TT=HH*0;
SS=HH*0;
nrec=0;
jj=578;
ii=1187;
CNC=HH*0;
for yr=yr1:yr2
  fmat = sprintf('%sarc04_%3.3i_precond_wint%i.mat',pthmat,expt,yr);
  fprintf('Loading %s\n',fmat);
  load(fmat);
  
  cnc=cnc+1;
  CONV=CONV+MXD.Conv_stat;
%  if MXD.Conv_stat(jj,ii)>0;
%    keyboard; 
%  end
  IT=find(~isnan(MXD.Tconv));
  TT(IT)=TT(IT)+MXD.Tconv(IT);
  SS(IT)=SS(IT)+MXD.Sconv(IT);
  CNC(IT)=CNC(IT)+1;
  nrec=nrec+MXD.Nrecords;
  
end

Prb=CONV./nrec;
Prb(Prb<0)=0;
Prb(Prb>1)=1;
IT=find(CNC>0);
TT(IT)=TT(IT)./CNC(IT);
SS(IT)=SS(IT)./CNC(IT);
TT(CNC==0)=nan;
SS(CNC==0)=nan;

c1=0;
c2=1;
%CMP = colormap_YOR(nint,c1,c2);
%cmp = CMP.colormap;
%for k=1:5
%  cmp(k,:)=[1 1 1];
%end
%cmp=smooth_colormap(cmp,15);
%cmp(1,:)=[1 1 1];
%cnt = CMP.intervals;
%nint=length(cmp);

xlim1=500;
xlim2=3200;
ylim1=800;
ylim2=3900;

% Probability
A.nf=1;
A.fld=Prb;
A.fldnm='prob';
A.xl1=xlim1;
A.xl2=xlim2;
A.yl1=ylim1;
A.yl2=ylim2;
A.c1c2=[0,1];
A.HH=HH;
A.LON=LON;
A.LAT=LAT;
A.title=sprintf('Cascading: %s 0.04-%3.3i %i-%i',...
  A.fldnm,expt,yr1,yr2);

sub_plot_precond(A);
set(gcf,'Position',[1480 359 1032 976]);
bottom_text(btx,'pwd',1);

% T
A.nf=2;
A.fld=TT;
A.fldnm='temp';
A.c1c2=[-2 8];
A.title=sprintf('Cascading: %s 0.04-%3.3i %i-%i',...
  A.fldnm,expt,yr1,yr2);

sub_plot_precond(A);
set(gcf,'Position',[1480 359 1032 976]);
bottom_text(btx,'pwd',1);

% S
A.nf=3;
A.fld=SS;
A.fldnm='salin';
A.c1c2=[26 36];
A.title=sprintf('Cascading: %s 0.04-%3.3i %i-%i',...
  A.fldnm,expt,yr1,yr2);

sub_plot_precond(A);
set(gcf,'Position',[1480 359 1032 976]);
bottom_text(btx,'pwd',1);

