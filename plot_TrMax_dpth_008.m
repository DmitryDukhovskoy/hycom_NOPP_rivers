% Plot depth of the maximum Tracer concentr
% or depth where integrated Tr.Mass exceeds 
% a threshold = % of the total depth-integrated
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

regn = 'ARCc0.08';
expt = 110;  

nTr   = 1;   % tracer to plot
s_fig = 0; % saves all active figures

TrMn=1e-2; % log Threshold value (conc), skip points<TrMn
rg = 9806;
fprintf('Tracer #: %i, Threhold value: %8.5d\n',nTr,TrMn);


rg = 9806;


%YRPLT=[2006,172];
YRPLT=[2016,365];
nyy = size(YRPLT,1);

np=size(YRPLT,1);

% Experiments:
pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc0.08/%3.3i/data_mat/',expt);


%ftopo = sprintf('%s/depth_%s_17DD.nc',pthtopo,regn); % 
ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2
MVOL = [];
[mm,nn] = size(HH);

for ip=1:nyy;
  yr   = YRPLT(ip,1);
  iday = YRPLT(ip,2);
  dnmb = datenum(yr,1,1)+iday-1;
  DV   = datevec(dnmb);
  
  fmat = sprintf('%s%s_%3.3i_VolIntgr_Tr%2.2i_%i%2.2i%2.2i.mat',...
		 pthmat,regn,expt,nTr,DV(1:3));
  
  
  %pthbin = sprintf('/nexsan/hycom/ARCc0.08_011/data/%i/',yr);  
  pthbin = sprintf('/nexsan/archive/%s_%3.3i/data/%i/',regn,expt,yr);
  fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
  finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);

  fprintf('\n:::: Analyzing %i/%2.2i/%2.2i   ::::\n\n',DV(1:3));

  fprintf('Computing vol-integrated Tracer mass ...\n');

  [F,n,m,l] = read_hycom(fina,finb,'tracer','r_tracer',nTr);
  F(F>1e6)=nan;
  Tr=F;

  [F,n,m,l] = read_hycom(fina,finb,'thknss');
  F=F./rg;
  F(F>1e10)=0;
  F(F<1e-2)=0;
  dH=F; 

  A = squeeze(Tr(1,:,:));
  IN = find(A>TrMn);
  
% Depth-integrate:
  Mtot = 0;
  Mtr = zeros(m,n); % mass per 1 m2, kg/m2
  for k=1:l
    fprintf('Integrating over vert. layers: k=%i\n',k);
    tr = squeeze(Tr(k,:,:)); % kg/m3
    dz = squeeze(dH(k,:,:));
    Mtr(IN) = Mtr(IN)+tr(IN).*dz(IN); % kg/m2 - normalized mass by area
  end
%keyboard
% Integrate to find depth where depth-intgr mass <= Cmax*Total-depth mass 
  Cmax = 0.8;
  cM = Mtr*0;
  Ib = IN;
  ln1 = length(IN);
  DTr = Mtr*0; % tracer depth
  for k=1:l
    tr = squeeze(Tr(k,:,:)); % kg/m3
    dz = squeeze(dH(k,:,:));
    cM(Ib) = cM(Ib)+tr(Ib).*dz(Ib); % kg/m2 - normalized mass by area
    rt = cM./Mtr;
    Ib = find(rt<=Cmax);
    ln2 = length(Ib);
    DTr(Ib) = DTr(Ib)+dz(Ib);
    rd = 1-ln2/ln1;
    fprintf('k=%2.2i: Finding max D tr>%3.2f, %6.4f%% pnts dropped \n',k,Cmax,rd*100);
  end
  
  DTr = - DTr;
  DTr(A<=TrMn)=nan;
  DTr(HH>-100)=nan;
  
%  xlim1 = 0;
%  xlim2 = nn;
%  ylim1 = 0;
%  ylim2 = mm;
% Greenland
  xlim1 = 380;
  xlim2 = 1250;
  ylim1 = 150;
  ylim2 = 1100;  
  nf = 1; % <0 - visible off
  ifx=max(strfind(fina,'/'));
  stl=sprintf('%s Depth of Mass=%4.2f*DpthIntgrMass, Tr=%i',...
		fina(ifx+1:end),Cmax,nTr);
  sub_plot_tracers(DTr,nf,HH,xlim1,...
		     xlim2,ylim1,ylim2,...
		     LON,LAT,stl,nTr,...
		     'c1',-1500,'c2',0,'cmp',4);


  txtb='plot_TrMax_dpth_008.m';
  bottom_text(txtb,'pwd',1,'position',[0.02 0.05 0.8 0.1]);

  if s_fig>0
    fgnm=sprintf('%sDpthMaxTrcr_tr%2.2i_%i%2.2i%2.2i',pthfig,nTr,DV(1:3));
    fprintf('Saving %s\n',fgnm);
    print('-dpng','-r250',fgnm);
  end
end
  




