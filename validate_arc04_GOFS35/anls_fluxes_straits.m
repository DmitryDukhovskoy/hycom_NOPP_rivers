% Summarize statistics for vol, heat, FW fluxes
% via main straits from the simulations
% extracted in extr_TSVdaily_straits04.m
%
% If >1 expt specifed - time ser. are combines into 1 filling gaps
%
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/Wavelet;
startup;

close all
clear

EXPT   = sub_cice_experiments;
iEX    = [9,6]; % use expts to have 1 full time series, first expt is the main uses 2nd to fill gaps
nmex1 = EXPT(iEX(1)).cice_opt;
nmex2 = EXPT(iEX(2)).cice_opt;


pthout = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.04/data_straits/';
%pthmat = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/%3.3i/strait_fluxes/',expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';


Tref1 = -1.8; % Ref T to calc. H flux
Tref2 = 0;    % Ref T to calc. H flux
Sref1 = 34.8;
Sref2 = 34.9;

f_zgrd=0;  % =1 - calculate fluxes from z-grid interpolated U,T,S - less accurate
           % mostly for comparison and validation

btx='anls_fluxes_straits.m';

ftopo = sprintf('%sdepth_ARCc0.04_17DD.nc',pthtopo); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);

% Fram Section is close to Moorings ~79N
%SCT = sub_define_sections04(HH,LON,LAT);
SCT = sub_define_AO_NA_sections04(HH,LON,LAT);
nsct = length(SCT);


CLR=[0 0.4 0.8; ...
     0.8 0.4 0; ...
     1 0.2 0; ...
     0 1 0; ...
     0.8 0 0.6; ...
     0.8 1 0; ...
     0.3 0.7 1; ...
     1 0.7 0.5; ...
     0.9 0.3 0.6; ...
     0.2 0.7 0.4; ...
     0.1 0.4 0.3; ...
     0.8 0.3 0.9; ...
     1 0.4 0.6; ...
     0.45 0.2 0.85];


f_map=0;
if f_map==1
  fprintf('Drawing map with segments\n');
  fn=1;
%  sub_plot_Greenl_contour(HH,LON,LAT,fn,GC);
  figure(1); clf;
  contour(HH,[0 0],'k');
  hold on;
  contour(HH,[-5000:500:-100],'Color',[0.9 0.9 0.9]);

  for ip=1:nsct
%    IJ=SCT(ip).IJ;
%    plot([IJ(1,1) IJ(2,1)],[IJ(1,2) IJ(2,2)],...
%        'Linewidth',2.5,'Color',[1. 0.6 0]);
   clr=CLR(ip,:);
    IIs=SCT(ip).I;
    JJs=SCT(ip).J;
    plot(IIs,JJs,'-',...
         'Linewidth',2.5,'Color',clr);
    Ip=SCT(ip).IJPR(1);
    Jp=SCT(ip).IJPR(2);

%    plot(Ip,Jp,'.','Markersize',14,'Color',clr);
  end

  axis('equal');
  set(gca,'xlim',[600 2800],...
          'ylim',[200 2800]);

  bottom_text(btx,'pwd',1);

  keyboard
end




FLX=[];
FLX(1).TM=[];
FLX(1).Vol=[];
FLX(1).FW1=[];
FLX(1).FW2=[];
YR=2017:2017
ix1 = iEX(1);
ix2 = iEX(2);

FLX = sub_combine2flx(ix1,ix2,EXPT,YR);

fprintf(' --------------------------------- \n\n');
nsct = length(FLX);
fplt = 1;
for ik=1:nsct
  nm=FLX(ik).Name;
  Vol=FLX(ik).Vol*1e-6;
%  FW1=FLX(ik).FW1*1e-3;
%  FW2=FLX(ik).FW2*1e-3;
  N=length(Vol);

  mvol=nanmean(Vol);
  p1=prctile(Vol,10);
  p9=prctile(Vol,90);
  verr=std(Vol)/sqrt(N); % mean error

%  ms1=nanmean(FW1);
%  p11=prctile(FW1,10);
%  p19=prctile(FW1,90);
%  s1err=std(FW1)/sqrt(N);
%
%  ms2=nanmean(FW2);
%  p21=prctile(FW2,10);
%  p29=prctile(FW2,90);
%  s2err=std(FW2)/sqrt(N);

  sv2km = 1e6*3600*24*365*1e-9; % Sv -> km3/yr
  msv2km = 1e3*3600*24*365*1e-9; % mSv -> km3/yr

  sinfo = sprintf('%s: VolFlx=%5.2f+/-%5.2f Sv, %6.1f+/-%6.1f km3/y3',...
    nm, mvol, verr, mvol*sv2km, verr*sv2km);

  if fplt==1
    figure(ik); clf;
    set(gcf,'Position',[1538         894        1001         424]);
    axes('Position',[0.09 0.4 0.84 0.4]);
    hold on;

    d1 = datenum(YR,1,1);
    tm1 = FLX(ik).TM1;
    tm2 = FLX(ik).TM2;
    tmo1 = days2months(tm1);
    tmo2 = days2months(tm2);	
    vol1 = FLX(ik).Vol1*1.e-6;
    vol2 = FLX(ik).Vol2*1.e-6;

    plot(tmo1,vol1,'-','Color',[1 0.3 0],'Linewidth',2.5);
    plot(tmo2,vol2,'-','Color',[0. 0.5 0.9],'Linewidth',2.5);
    
    set(gca,'Tickdir','out',...
            'xlim',[1 ceil(max([tmo1(end),tmo2(end)]))],...
            'xtick',[0:ceil(max([tmo1(end),tmo2(end)]))],...
            'Fontsize',14,...
            'xgrid','on',...
            'ygrid','on');

    lgd = legend(nmex1,nmex2);
    set(lgd,'Position',[0.7 0.18 0.2 0.11])
    stl = sprintf('VolFlx, %i, %s',YR,nm);
    title(stl);

    axes('Position',[0.05,0.2,0.4,0.05]);
    text(0.2,0.5,sinfo,'Fontsize',12);
    set(gca,'ylim',[0.3 0.6],'visible','off');


    drawnow
  end

%  fprintf('%s: Vol=%5.2f+/-%5.2f Sv, %6.1f+/-%6.1f km3/y3\n',...
%    nm, mvol, verr, mvol*sv2km, verr*sv2km);
  fprintf('%s\n',sinfo);

  fprintf(' --------------------------------- \n\n');

end






