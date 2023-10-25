% Plot volume integrated mass of the tracer
% extracted in vol_intgr_trcr008_monthly.m
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

regn = 'ARCc0.08';
expt = 110;  
YR1 = 1993;
YR2 = 2015;

nTr   = 1   % tracer to plot

% Experiments:
pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc0.08/%3.3i/data_mat/',expt);
btx = 'plot_vol_intgr_trcr008_monthly.m';

MVOL = [];
ZM = [];
TM = [];
for iyr=YR1:YR2
  yr   = iyr;
  fmat = sprintf('%s%3.3i_VrtLrMass_mo_Tr%2.2i_%i.mat',...
		 pthmat,expt,nTr,iyr);
  fprintf('Loading %s\n',fmat);
  load(fmat);
  
  tm = TRI.TM;
  tm = tm(:)';
  TM = [TM,tm];
% Mass tracers by layers:
  dzm  = TRI.Mean_LThkn;
  zzm  = -cumsum(dzm);
  Mv   = TRI.Vert_LrMass_kg*1e-12; % kg->tonn (1e-3)->GT (1e-9) in each layer
% Layers have different thickness, need to normalize to compare
% thus use Mass/1 m of depth
  Mvm  = Mv./abs(dzm);
  ZM   = [ZM,zzm];
  MVOL = [MVOL,Mvm]; % GT (1e-9)
%  keyboard
end; % time loop

DV  = datevec(TM);
yrs = [DV(1,1):1/12:DV(end,1)+0.99];
[YRS,dmm]=meshgrid(yrs,ZM(:,1));

%nint=200;
%c1=0;
%c2=6;
%CMP = create_colormap5(nint,c1,c2);
%CMP = create_colormap4(nint,c1,c2);
%cnt = CMP.intervals;
%cmp = CMP.colormap;

nint = 320;
c1 = 0;
c2 = 8;
%CMP = create_colormap2_3(nint,c1,c2);
CMP = colormap_sclr2(nint,c1,c2);
cmp = CMP.colormap;
for ik=1:20
  cmp(ik,:) = [1 1 1];
end
cmp = smooth_colormap(cmp,20);
cmp = smooth_colormap(cmp,20);
cmp(1,:) = [1 1 1];
cnt = CMP.intervals;


figure(1); clf;
pcolor(YRS,ZM,MVOL); shading interp;
caxis([c1 c2]);
colormap(cmp);

hc=colorbar;

bottom_text(btx,'pwd',1);



