% Plot time series of the annual dlt(S)
% estimated from 2 epxeriments: with and without Gr.Flux
% S is averaged for several v. layers
% and extracted in mnthly_arc08_layers_S.m
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

s_mat = 2; 
YR1  = 1993;
YR2  = 2016;
imo = 13; % 13 - yearly average

LRS = load('LRS.dat');
nlrs= length(LRS)-1; % skip whole depth  
       
regn = 'ARCc0.08';
expt = 110; % no Greenland runoff  
%expt = 112;  % Greenland runoff

s_fig = 0;

% Experiments:
% 102 - test simulation, HYCOM GLBb0.08 nest: 1 file kept constant (1993,1,1)

pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
%pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc0.04/%3.3i/data_mat/',expt);

fmat = sprintf('%sarc08_%3.3i_dS_tser_ATltRegions.mat',pthmat,expt);

%ftopo = sprintf('%s/depth_%s_17DD.nc',pthtopo,regn); % 
ftopo = sprintf('%s/depth_%s_09.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

%[DX,DY]=sub_dx_dy(LON,LAT);
%Acell=DX.*DY; % Grid cell area, m2

BX = sub_define_boxes(HH,LON,LAT,0);

[XX,YY] = meshgrid((1:nn),(1:mm));
for ib=1:5
  iBG = BX(ib).IJ;
  INP = inpolygon(XX,YY,iBG(:,1),iBG(:,2));
  IN = find(INP==1);
  BX(ib).IN = IN;
end

% Butterworth filter
% output freq. 1/12 mo^-1
% 12 mo cutoff = 12/12=1 -> in Matlab Wn=1/6
% 
Wn = 1/6; % cutoff freq 6 mo: 2/6, 1yr=1/6
[Bf,Af] = butter(9,Wn,'low');


if s_mat==1
  cc=0;
  for iyr=YR1:YR2
    YR=iyr;
    for mo=1:12
      cc=cc+1;
      expt1=110;
      pthmat1 = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt1);
      pthm1 = sprintf('/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_mat%3.3i/',expt1);
      fmat1 = sprintf('%sarc08_%3.3i_mnthS_lrs_%4.4i%2.2i.mat',pthm1,expt1,YR,mo);
      fprintf('Loading %s\n',fmat1);
      load(fmat1);
      S1 = meanS;

      expt2=112;
      pthmat2  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt2);
      pthm2 = sprintf('/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_mat%3.3i/',expt2);
      fmat2 = sprintf('%sarc08_%3.3i_mnthS_lrs_%4.4i%2.2i.mat',pthm2,expt2,YR,mo);
      fprintf('Loading %s\n',fmat2);
      load(fmat2);
      S2 = meanS;

      for ilv=1:4
	ss1=S1(ilv).Savrg;
	ss2=S2(ilv).Savrg;
	dS(ilv,:,:) = ss2-ss1; 
      end

      for ib=1:5
	nm=BX(ib).Name;
	IN=BX(ib).IN;
	for ilv=1:4
	  dss=squeeze(dS(ilv,:,:));
	  dmm=nanmean(dss(IN));
	  DltS(ib,ilv,cc)=dmm;
	end
      end

    end
  end
  fprintf('Saving %s\n',fmat);
  save(fmat,'DltS');
else
  fprintf('Loading %s\n',fmat);
  load(fmat);
end  

%keyboard

%CLR = [0.9 .9 0.9; ...
%       0.6 0.6 0.6; ...
%       0.45 0.45 0.45; ...
%       0.1 0.1 0.1];
CLR = [0 0.4 0.6; ...
       0.7 0.3 0; ...
       0.8 0.7 0.; ...
       0. 0.8 0.4];

yrs=[1993:1/12:2016.99];
for ib=1:5
  nm=BX(ib).Name;
  figure(ib); clf;
  axes('Position',[0.09 0.45 0.85 0.45]);
  hold on;
  for ilv=1:4
    dS=squeeze(DltS(ib,ilv,:));
    dmm = filtfilt(Bf,Af,dS);
    clr = CLR(ilv,:);
    plot(yrs,dmm,'Linewidth',1.6,'Color',clr);
  end
  
  set(gca,'tickdir','out',...
	  'xlim',[1993 2017],...
	  'xtick',[1993:2018],...
	  'xgrid','on',...
	  'ygrid','on',...
	  'Fontsize',12);
	  

  hll=legend('0-50','50-150','150-300','300-500',...
	 'Location','SouthOutside');
  set(hll,'Position',[0.52 0.26 0.1 0.1],'Fontsize',12);

  stl=sprintf('dS from "No/With GrRunoff" expt, %s',nm);
  title(stl);
  
  btx='dltS_tser_LrAvrg_GreenldExp_NatlRegions.m';
  bottom_text(btx,'pwd',1,'position',[0.02 0.3 0.4 0.1]);
  
end

  




