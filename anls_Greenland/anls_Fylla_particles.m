% Analyze lagr. particles pathways
% to determine Greenland FW pathways
% for calculating propagation speed
% see Fylla_particles.m
%
% 2 strategies: use daily U,V - lots of eddies
% use annual mean U,V - straight pathways
%  initialize with different meanUV fields
% to get different pathways

addpath /Net/Movies0/ddmitry/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /home/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers;
%addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
%addpath /Net/Movies0/ddmitry/MyMatlab/seawater
startup

close all
clear

ufld = 'mean'; % mean - annual mean fields, daily - daily U
regn = 'ARCc0.08';
expt = 110;  
YR0 = 2005; % for mean UV - year when particles initialized
YR1 = YR0;
YR2 = 2010;

pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/110/data_Gr_prt/';

% ------------------------
% TOPO
% ------------------------
ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

hmsk=HH;
hmsk(HH<0)=nan;

IP=[];
JP=[];
TM=[];
for ir=YR1:YR2
  YR=ir;
  switch(ufld),
    case('daily');
     fmat = sprintf('%sFylla_particles_%i.mat',pthmat,YR); % with daily fields
    case('mean');
  fmat = sprintf('%sFylla_particles_meanUyr%i_%i.mat',pthmat,YR0,YR); % with daily fields
  end
  
  fprintf('\n   Loading saved %s\n\n',fmat);
  load(fmat);
  TR=PRTCL.TRACK;

  nt = length(TR);
  % If there are several time steps saved, plot tracks"
  TR(1).TM=datenum(1993,1,1);
  clear tm
  for it=1:nt
    tm=TR(it).TM;
    if ~isempty(TM) & tm==TM(end), continue; end;    
    TM=[TM;tm];
    dmm=TR(it).I;
    IP=[IP,dmm];
    dmm=TR(it).J;
    JP=[JP,dmm];
  end

end; % years

figure(1); clf;
contour(HH,[0 0],'k','Linewidth',1.6);
hold on;
contour(HH,[-200 -200],'Color',[0.8 0.8 0.8]);
contour(HH,[-5000:1000:-100],'Color',[0.6 0.6 0.6]);

npp=size(IP,1);
for ipt=1:npp
  plot(IP(ipt,:),JP(ipt,:),'-');
end

axis('equal');
set(gca,'xlim',[350 1250],...
	'ylim',[1 1100]);

stl=sprintf('Fylla Bank, traj, U %s, %i-%i',ufld, YR1,YR2);