% Analyze S tseries 
% extracted in extr_S_tser_vertlrs_pnt.m
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

YR1=2006;
YR2=2006;
regn = 'ARCc0.08';
%expt = 110; % no Greenland runoff  
expt = 112;  % Greenland runoff

fprintf('expt %3.3i\n\n',expt);


s_mat = 1; % =2 - load and start from last saved
%pfld  = 'salin';

s_fig = 0;
rg = 9806;
fprintf('T&S, %i-%i\n',YR1,YR2);


pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);

btx='plot_Stser.m';

ftopo = sprintf('%s/depth_%s_09.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

% Combine time series from several years:
TM=[];
SS=[];
TT=[];
ZZ=[];
IJp=[];
TM=[];
for yr=YR1:YR2
  fmat = sprintf('%sarc08_%3.3i_TS_timeser_%i.mat',pthmat,expt,yr);
  fprintf('Loading %s\n',fmat);
  load(fmat);
  
  if isempty(IJp)
    IJp=TS.IJ_points;
  end
  
  Z=TS.Z;
  S=TS.S;
  T=TS.T;
  tm=TS.TM;
  ZZ=[ZZ;Z];
  SS=[SS;S];
  TT=[TT;T];
  TM=[TM;tm(:)];
end



% Points, arc08 indices:
IJp=[490         596
         431         449
         509         401
         657         176
         712         470
         944         418
         967         611
        1017         790
        1079         546
        1151         690];
np=length(IJp);

% Plot
f_chck=0;
if f_chck==1
  figure(10); clf;
  contour(HH,[0 0],'k');
  hold on;
  contour(HH,[-10000:1000:-900],'color',[0.7 0.7 0.7]);
  axis('equal');
  
  for ik=1:np
    i0=IJp(ik,1);
    j0=IJp(ik,2);
    plot(i0,j0,'r*');
    text(i0+5,j0,sprintf('%2.2i',ik),'Fontsize',12);
  end
  title('T/S time series locations');
  
  bottom_text(btx,'pwd',1);
end


DV=datevec(TM);
Td=(TM-TM(1))/365;
nlrs=size(SS,2);
ntm=length(TM);
for ip=1:np;
  Z=squeeze(ZZ(ip,:,:));
  S=squeeze(SS(ip,:,:));
  
  
  dZ=zeros(nlrs,1);
  Sav=zeros(nlrs,1);
  ss0=squeeze(SS(ip,1,:));
% subtract seasonal mean:
  Smnth=[];
  for im=1:12
    Im=find(DV(:,2)==im);
    a=ss0(Im);
    ss0(Im)=ss0(Im)-mean(a);
  end
  dltS(1,:)=ss0;
  
 
  for k=1:nlrs
    dz=squeeze(abs(ZZ(ip,k+1,:)-ZZ(ip,k,:)));
    ss=squeeze(SS(ip,k,:));
    Sav(k,:)=nansum(ss.*dz);
    dZ(k,:)=nansum(dz);
    ss1(k,:)=prctile(ss,10);
    ss2(k,:)=prctile(ss,90);
%
% Correlation with surface
% subtract seasonal mean:
    for im=1:12
      Im=find(DV(:,2)==im);
      a=ss(Im);
      ss(Im)=ss(Im)-mean(a);
      Smnth(k,im)=mean(a);
    end
    dltS(k,:)=ss;
    nlg=10;
    CC=xcorr(ss0,ss,'coeff',nlg);
    Crr(k)=CC(nlg+1);
  end
  Sav=Sav./dZ;
  dZ=dZ./ntm; % average layer thickness
  
  
% Plot profiles:  
  Sav(dZ==0)=nan;
  ss1(dZ==0)=nan;
  ss2(dZ==0)=nan;
  
  Zav=-cumsum(dZ);
  xl1=min(min(Smnth));
  xl2=max(max(Smnth(1:26,:)));
  yl1=-200;
  
  
  f_plt1=1;  % plot individual locations on each fig.
  if f_plt1==1
  
  figure(ip); clf;
% 
% August profile:
    axes('Position',[0.09 0.4 0.23 0.5]);
    plot(Smnth(:,8),Zav,'Linewidth',2,'Color',[0.9 0.3 0]); 
    hold on;
    plot(Smnth(:,3),Zav,'Linewidth',2,'Color',[0 0.5 0.8]); 
    set(gca,'tickdir','out',...
	    'xlim',[xl1 xl2],...
	    'ylim',[yl1 0],...
	    'xgrid','on',...
	    'ygrid','on',...
	    'ytick',[-400:25:0],...
	    'Fontsize',14);
    stl=sprintf('0.08-112, Pnt %i, Aug & March S, %i-%i',ip,YR1,YR2);
    title(stl);


    axes('Position',[0.42 0.4 0.23 0.5]);
    plot(Sav,Zav,'Color',[0 0.3 0.7],'Linewidth',2.2); 
    hold on;
    plot(ss1,Zav,'-','Linewidth',2,'Color',[0. 0.6 1]);
    plot(ss2,Zav,'-','Linewidth',2,'Color',[0. 0.6 1]);
    set(gca,'tickdir','out',...
	    'ylim',[yl1 0],...
	    'xgrid','on',...
	    'ygrid','on',...
	    'ytick',[-400:25:0],...
	    'Fontsize',14);
    stl=sprintf('mean S',ip,YR1,YR2);
    title(stl);

  % X-correlation
    axes('Position',[0.75 0.4 0.23 0.5]);
    plot(Crr,Zav,'Linewidth',2); 
    set(gca,'tickdir','out',...
	    'xlim',[min(Crr)-0.02 1.02],...
	    'ylim',[yl1 0],...
	    'xgrid','on',...
	    'ygrid','on',...
	    'ytick',[-200:25:0],...
	    'xtick',[-1:0.2:1],...
	    'Fontsize',14);
    stl=sprintf('XCorr Surf-S(z)');
    title(stl);

    bottom_text(btx,'pwd',1);

  end
  
  
% Plot all profiles on 1 figure  
  f_all=0;
  if f_all==1
    if ~exist('CLR','var')
      CLR=rand(np,3);
    end

    figure(np+1);
    if ip==1
      clf;
      axes('Position',[0.12 0.4 0.25 0.5]);
      hold on
    end
    clr=CLR(ip,:);
    plot(Smnth(:,8),Zav,'Linewidth',2,'Color',clr); 

    set(gca,'tickdir','out',...
	    'ylim',[yl1 0],...
	    'xgrid','on',...
	    'ygrid','on',...
	    'ytick',[-400:25:0],...
	    'Fontsize',14);
    title('August S');
    if ip==np
      bottom_text(btx,'pwd',1);
    end

    figure(np+2);
    if ip==1
      clf;
      axes('Position',[0.12 0.4 0.25 0.5]);
      hold on
    end
    plot(Smnth(:,3),Zav,'Linewidth',2,'Color',clr); 

    set(gca,'tickdir','out',...
	    'ylim',[yl1 0],...
	    'xgrid','on',...
	    'ygrid','on',...
	    'ytick',[-400:25:0],...
	    'Fontsize',14);
    title('March S');
    if ip==np
      bottom_text(btx,'pwd',1);
    end

  % Correlation  
    figure(np+3);
    if ip==1
      clf;
      axes('Position',[0.12 0.4 0.25 0.5]);
      hold on
    end

    plot(Crr,Zav,'Linewidth',2,'Color',clr); 
    set(gca,'tickdir','out',...
	    'ylim',[yl1 0],...
	    'xgrid','on',...
	    'ygrid','on',...
	    'ytick',[-400:25:0],...
	    'Fontsize',14);
    title('Correlation S0-S(z)');
    if ip==np
      bottom_text(btx,'pwd',1);
    end

    figure(np+4);
    if ip==1
      clf;
      axes('Position',[0.12 0.4 0.25 0.5]);
      hold on
    end

    yy0=np-ip;
    plot([0 1],[yy0 yy0],'-','Linewidth',2,'Color',clr); 
    text(1.2,yy0, sprintf('Pnt %i',ip),'Fontsize',14);

    if ip==np
      set(gca,'xlim',[0 4],...
	      'visible','off',...
	      'Fontsize',14);
    end
  
  end
  
end





