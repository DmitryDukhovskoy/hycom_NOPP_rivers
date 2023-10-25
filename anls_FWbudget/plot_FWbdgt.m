% Caculate FW content wrt to Sref=34.8
% full depth and to isohaline = Sref (Haine et al, 2015)
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

expt = 110;
s_fig  = 1;

rg=9806;  % convert pressure to depth, m
Sref=34.8;

regn = 'ARCc0.08';
btx = 'plot_FWbdgt.m';

%pthbin  = '/nexsan/GLBb0.08/';
pthmat  = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/data_mat/';
pthfig  = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/fig_fwc/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';

fmat    = sprintf('%sfwc_hycom_arc08_%3.3i.mat',pthmat,expt);

ftopo = sprintf('%sdepth_ARCc0.08_11.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(HH);


fprintf('Loading %s\n',fmat);
load(fmat);

nb=length(FWBX);
for ik=nb-1:nb-1rc08_110_fwc_ArcticOcean.png
  nm   = FWBX(ik).Name;
  TM   = FWBX(ik).TM;
  fwc1 = FWBX(ik).Fwc1_m; % depth-integrated FWC
  fwc2 = FWBX(ik).Fwc2_m; % integrated from Sref depth
  
  DV=datevec(TM); % monthly data
  nc=length(DV);
  nyr=nc/12;
  AA1 = reshape(fwc1,[12,nyr]);
  am1 = mean(AA1);
  AA2 = reshape(fwc2,[12,nyr]);
  am2 = mean(AA2);
  TT=[0:nc-1]'/12+DV(1,1);
  TY=[DV(1,1):DV(end,1)];

%  figure('Visible','off'); clf;
  figure(1); clf;
%  set(gcf,'Visible','off');
  axes('Position',[0.1 0.58 0.8 0.35]);
%  bb=bar(TY,am1-am1(1));
  bb=bar(TY,am1,0.92);
  set(bb,'Facecolor',[0. 0.4 0.7],'EdgeColor','none');
  sll=sprintf('%s, Dpth-Intgr, FWC Srf=%4.1f',...
	      nm,FWBX(ik).Sref);
  title(sll);
  set(gca,'tickdir','out',...
	  'xlim',[DV(1,1)-0.5 DV(end,1)+0.5],...
	  'xgrid','on',...
	  'xtick',[1991:2017]);
  
  axes('Position',[0.1 0.09 0.8 0.35]);
%  bb=bar(TY,am2-am2(1));
%  set(bb,'Facecolor',[0. 0.4 0.7],'EdgeColor','none');
%  sll=sprintf('%s, d(FWC) wrt 1993=%4.1f m, Intgr z=Srf=%4.1f',...
%	      nm,am2(1),FWBX(ik).Sref);
%  title(sll);
  bb=bar(TY,am2,0.92);
  set(bb,'Facecolor',[0. 0.4 0.7],'EdgeColor','none');
  sll=sprintf('%s, FWC(m) Intgr z=Srf=%4.1f',...
	      nm,FWBX(ik).Sref);
  title(sll);
  set(gca,'tickdir','out',...
	  'xlim',[DV(1,1)-0.5 DV(end,1)+0.5],...
	  'xgrid','on',...
	  'xtick',[1991:2017]);
  
  bottom_text(btx,'pwd',1);
  if s_fig==1
    fgnm = sprintf('%sarc08_110_fwc_%s',pthfig,nm);
    fprintf('Saving %s\n',fgnm);
    print('-dpng','-r250',fgnm);
  end
end  
  
  
  
  
