% Plot fluxes
% for POP sections on the Greenland shelf
% Fluxes calculated in greenl_fluxes_POPgates.m
%
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/Wavelet;
startup;

close all
clear


expt=110;
TV=11;
YR1=2005;
YR2=2008;


pthmat = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_theresa/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';

btx='plot_greenl_fluxes_POPgates.m';

fprintf('Plotting arc08-%3.3i fluxes Greenland Shelf gates %i-%i\n',...
	expt,YR1,YR2);


ftopo = sprintf('%sdepth_ARCc0.08_%2.2i.nc',pthtopo,TV); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);

GC = sub_greenl_isobath(HH,LON,LAT);
Ig=GC.cntr_Iindx;
Jg=GC.cntr_Jindx;

%
% Combine time series:
FLX=struct;
cc=0;
for YR=YR1:YR2;
  cc=cc+1;
  fmatout=sprintf('%shycom008_%3.3i_Greenl_flx_POPgates_%4.4i.mat',...
	      pthmat,expt,YR);
  fprintf('Load %s\n',fmatout);
  load(fmatout);
  if cc==1
    for ik=1:6
      FLX(ik).Name=VHFLX(ik).Name;
      FLX(ik).T1=VHFLX(ik).Tref1;
      FLX(ik).T2=VHFLX(ik).Tref2;
      FLX(ik).Vol=VHFLX(ik).VolFlx_m3s;
      FLX(ik).Hf1=VHFLX(ik).HFlx_T1_W;
      FLX(ik).Hf2=VHFLX(ik).HFlx_T2_W;
%      FLX(ik).Time=VHFLX(ik).Time;
      I=length(VHFLX(ik).VolFlx_m3s);
      yday=YR+[0:I-1]./I;
      FLX(ik).Time=yday;
      FLX(ik).VolGr=VHFLX(ik).VolFlxGrSh_m3s;
      FLX(ik).Hf1Gr=VHFLX(ik).HFlxGrSh_T1_W;
      FLX(ik).Hf2Gr=VHFLX(ik).HFlxGrSh_T2_W;
    end
    
  else
    for ik=1:6
      dmm=[FLX(ik).Vol,VHFLX(ik).VolFlx_m3s];
      FLX(ik).Vol=dmm;
      dmm=[FLX(ik).Hf1,VHFLX(ik).HFlx_T1_W];
      FLX(ik).Hf1=dmm;
      dmm=[FLX(ik).Hf2,VHFLX(ik).HFlx_T2_W];
      FLX(ik).Hf2=dmm;
%      dmm=[FLX(ik).Time,VHFLX(ik).Time];
%      FLX(ik).Time=dmm;    
      I=length(VHFLX(ik).VolFlx_m3s);
      yday=YR+[0:I-1]./I;
      dmm=[FLX(ik).Time,yday];
      FLX(ik).Time=dmm;
      dmm=[FLX(ik).VolGr,VHFLX(ik).VolFlxGrSh_m3s];
      FLX(ik).VolGr=dmm;
      dmm=[FLX(ik).Hf1Gr,VHFLX(ik).HFlxGrSh_T1_W];
      FLX(ik).Hf1Gr=dmm;
      dmm=[FLX(ik).Hf2Gr,VHFLX(ik).HFlxGrSh_T2_W];
      FLX(ik).Hf2Gr=dmm;
    end
  end
    
end


% Itegrated over all section fluxes:
for ik=1:6
  nm=FLX(ik).Name;
  TMd=FLX(ik).Time;
  T1=VHFLX(1).Tref1;
  T2=VHFLX(1).Tref2;
  
  VF=FLX(ik).Vol*1e-6;  % Sv
  mVF=mean(VF);
  pVF10=prctile(VF,10);
  pVF90=prctile(VF,90);
  
  HF1=FLX(ik).Hf1*1e-12;
  mHF1=mean(HF1);
  pHF110=prctile(HF1,10);
  pHF190=prctile(HF1,90);
  
  HF2=FLX(ik).Hf2*1e-12;
  mHF2=mean(HF2);
  pHF210=prctile(HF2,10);
  pHF290=prctile(HF2,90);

  
  figure(ik); clf;
  axes('Position',[0.1 0.7 0.8 0.22]);
  hold on
  plot(TMd,VF);
  plot([TMd(1) TMd(end)],[mVF mVF],'r--');
  plot([TMd(1) TMd(end)],[pVF10 pVF10],'--','Color',[0.6 0.6 0.6]);
  plot([TMd(1) TMd(end)],[pVF90 pVF90],'--','Color',[0.6 0.6 0.6]);
  
  set(gca,'xlim',[YR1 YR2+0.01],...
	  'ylim',[1.01*min(VF) 1.01*max(VF)],...
	  'xtick',[YR1:0.5:YR2],...
	  'xgrid','on',...
	  'ygrid','on');
  
  stl=sprintf('HYCOM 008-110, %s, VolFlux, %4.2f Sv %i-%i',nm,mVF,YR1,YR2);
  title(stl);
  
% Heat Flux 1  
  axes('Position',[0.1 0.4 0.8 0.22]);
  hold on
  plot(TMd,HF1);
  plot([TMd(1) TMd(end)],[mHF1 mHF1],'r--');
  plot([TMd(1) TMd(end)],[pHF110 pHF110],'--','Color',[0.6 0.6 0.6]);
  plot([TMd(1) TMd(end)],[pHF190 pHF190],'--','Color',[0.6 0.6 0.6]);
  
  set(gca,'xlim',[YR1 YR2+0.01],...
	  'ylim',[1.01*min(HF1) 1.01*max(HF1)],...
	  'xtick',[YR1:0.5:YR2],...
	  'xgrid','on',...
	  'ygrid','on');
  
  stl=sprintf('HYCOM 008-110, %s, HeatFlx, %5.2f TW, T0=%3.1f, %i-%i',...
	      nm,mHF1,T1,YR1,YR2);
  title(stl);
  
% Heat Flux 2  
  axes('Position',[0.1 0.1 0.8 0.22]);
  hold on
  plot(TMd,HF2);
  plot([TMd(1) TMd(end)],[mHF2 mHF2],'r--');
  plot([TMd(1) TMd(end)],[pHF210 pHF210],'--','Color',[0.6 0.6 0.6]);
  plot([TMd(1) TMd(end)],[pHF290 pHF290],'--','Color',[0.6 0.6 0.6]);
  
  set(gca,'xlim',[YR1 YR2+0.01],...
	  'ylim',[1.01*min(HF2) 1.01*max(HF2)],...
	  'xtick',[YR1:0.5:YR2],...
	  'xgrid','on',...
	  'ygrid','on');
  
  stl=sprintf('HYCOM 008-110, %s, HeatFlx, %5.2f TW, T0=%3.1f, %i-%i',...
	      nm,mHF2,T2,YR1,YR2);
  title(stl);
  
  
  
end

bottom_text(btx,'pwd',1);


% Plot fluxes only on Gr. Section:
for ik=1:6
  nm=FLX(ik).Name;
  TMd=FLX(ik).Time;
  T1=VHFLX(1).Tref1;
  T2=VHFLX(1).Tref2;
  
  VF=FLX(ik).VolGr*1e-6;  % Sv
  mVF=mean(VF);
  pVF10=prctile(VF,10);
  pVF90=prctile(VF,90);
  
  HF1=FLX(ik).Hf1Gr*1e-12;
  mHF1=mean(HF1);
  pHF110=prctile(HF1,10);
  pHF190=prctile(HF1,90);
  
  HF2=FLX(ik).Hf2Gr*1e-12;
  mHF2=mean(HF2);
  pHF210=prctile(HF2,10);
  pHF290=prctile(HF2,90);

  
  figure(ik+10); clf;
  axes('Position',[0.1 0.7 0.8 0.22]);
  hold on
  plot(TMd,VF);
  plot([TMd(1) TMd(end)],[mVF mVF],'r--');
  plot([TMd(1) TMd(end)],[pVF10 pVF10],'--','Color',[0.6 0.6 0.6]);
  plot([TMd(1) TMd(end)],[pVF90 pVF90],'--','Color',[0.6 0.6 0.6]);
  
  set(gca,'xlim',[YR1 YR2+0.01],...
	  'ylim',[1.01*min(VF) 1.01*max(VF)],...
	  'xtick',[YR1:0.5:YR2],...
	  'xgrid','on',...
	  'ygrid','on');
  
  stl=sprintf('HYCOM 008-110, %s, Gr VFlx, %4.2f Sv %i-%i',nm,mVF,YR1,YR2);
  title(stl);
  
% Heat Flux 1  
  axes('Position',[0.1 0.4 0.8 0.22]);
  hold on
  plot(TMd,HF1);
  plot([TMd(1) TMd(end)],[mHF1 mHF1],'r--');
  plot([TMd(1) TMd(end)],[pHF110 pHF110],'--','Color',[0.6 0.6 0.6]);
  plot([TMd(1) TMd(end)],[pHF190 pHF190],'--','Color',[0.6 0.6 0.6]);
  
  set(gca,'xlim',[YR1 YR2+0.01],...
	  'ylim',[1.01*min(HF1) 1.01*max(HF1)],...
	  'xtick',[YR1:0.5:YR2],...
	  'xgrid','on',...
	  'ygrid','on');
  
  stl=sprintf('HYCOM 008-110, %s, Gr HFlx, %5.2f TW, T0=%3.1f, %i-%i',...
	      nm,mHF1,T1,YR1,YR2);
  title(stl);
  
% Heat Flux 2  
  axes('Position',[0.1 0.1 0.8 0.22]);
  hold on
  plot(TMd,HF2);
  plot([TMd(1) TMd(end)],[mHF2 mHF2],'r--');
  plot([TMd(1) TMd(end)],[pHF210 pHF210],'--','Color',[0.6 0.6 0.6]);
  plot([TMd(1) TMd(end)],[pHF290 pHF290],'--','Color',[0.6 0.6 0.6]);
  
  set(gca,'xlim',[YR1 YR2+0.01],...
	  'ylim',[1.01*min(HF2) 1.01*max(HF2)],...
	  'xtick',[YR1:0.5:YR2],...
	  'xgrid','on',...
	  'ygrid','on');
  
  stl=sprintf('HYCOM 008-110, %s, Gr HFlx, %5.2f TW, T0=%3.1f, %i-%i',...
	      nm,mHF2,T2,YR1,YR2);
  title(stl);
  
  
  
end

bottom_text(btx,'pwd',1);

