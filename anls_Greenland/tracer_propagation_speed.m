% To estimate propagation speed of FW/tracer
% extract Tracer conc along specified pathways
% see anls_UV/tracer_pathways_meanUV.m
% pathways - using spline interpolation
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

s_mat=2; % =0 - extract data, no save; =1 extract & save, =2 - load saved
s_fig=0;
rg=9806;  % convert pressure to depth, m
btx = 'tracer_propagation_speed.m';

nTr = 1; % what tracer
npth = 1; % pathway #

YRPLT=[];
cc=0;
iyr2=2010;
for iyr=1993:2016
  for im=1:12
    cc=cc+1;
    YRPLT(cc,1)=iyr;
    YRPLT(cc,2)=im;
  end
end

%YRPLT=[2008,50];
np=size(YRPLT,1);


regn = 'ARCc0.08';
expt = 110;
pthfig  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/fig_sections/',expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);

fmpts = sprintf('%strcr_pathway_%2.2i.mat',pthmat,npth);
fmout = sprintf('%strcr_propagation_NAtl_pthw%2.2i.mat',pthmat,npth);


fprintf('Plotting xsection for %s-%3.3i, Tracer %2.2i\n\n',regn,expt,nTr);

%ftopo = sprintf('%s/depth_%s_09.nc',pthtopo,regn); % 
ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

% Specify segments of the x-section
fprintf('Loading pathway points: %s\n',fmpts);
load(fmpts);

IJs = PTS.IJpathway;
IIs = IJs(:,1);
JJs = IJs(:,2);
nS  = length(IIs);
Dst = PTS.Dist;
INDs= PTS.Ipathway;

f_map=0;
if f_map>0
  figure(10); clf;
  axis('equal');
  hold on
  contour(HH,[0 0],'k');
  contour(HH,[-5000:1000:-100],'Color',[0.7 0.7 0.7]);
  plot(IIs,JJs,'b.-');
  for ik=1:11
    xx=(ik-1)*1000;
    D=abs(Dst-xx);
    i0=find(D==min(D),1);
    plot(IIs(i0),JJs(i0),'r.','Markersize',14);
    stl=sprintf('%i',xx);
    text(IIs(i0),JJs(i0),stl,'Fontsize',12);
  end
  
  bottom_text(btx,'pwd',1);
end

% Get tracers, from monthly mean
% extracted for depth layers
if s_mat==1
  TRC = [];
  for ip=1:np
    iyr=YRPLT(ip,1);
    imo=YRPLT(ip,2);

    fmat = sprintf('%strcr_dpthav_%i%2.2i.mat',pthmat,iyr,imo);
    fprintf('Loading %s\n',fmat);
    if exist(fmat,'file')
      load(fmat);
    else
      fprintf(' =========  MISSING %s\n',fmat);
      return
    end
    rr = [];

    for ilv=1:3
      fprintf('Reading: %i/%2.2i, Tracer %i, Lev %i\n',...
	    iyr,imo,nTr,ilv);

      Tr = squeeze(TRCR(ilv).TR(nTr,:,:));
      Tr(Tr<=1e-23)=nan;

      TRC(ilv,ip,:)=Tr(INDs);
    end
  end

  if s_mat==1
    fprintf('Saving %s\n',fmout);
    save(fmout,'TRC');
  end

else
  fprintf('Loading %s\n',fmout);
  load(fmout);
end

clear PF
xyr=[1993:1/12:2016.99];
for ilv=1:2
  aa=squeeze(TRC(ilv,:,:));
  aa(aa<1e-8)=nan;
  laa=log(aa);

% Find when conc =80% of the max at every point
  [ma,na]=size(aa);
  ny=ma/12;
%  prt=0.8;
  clear YRTR  % keep years when tracer conc reaches quasi-steady state
  for ik=1:na
    a=aa(:,ik);
    ay=reshape(a,12,ny);
    ava=nanmean(ay);  % annual mean for current location
    ava(ava==0)=nan;
    dav=diff(ava);
% find year when the rest of the years the trend change is ~0 trend
    alf=100;
    Eps=max(ava)/240;
    iR=[];
    ir1=min(find(~isnan(ava)));
%    if ir1==1, ir1=2; end;
    BB=zeros(ny,2)*nan;
    YF=zeros(ny,ny);
    for iir=ir1:ny-1
      bb=ava(iir:end)';
      nb=length(bb);
      XX=[ones(nb,1),[iir:ny]'];
      B=regress(bb,XX);
      BB(iir,1)=B(1);
      BB(iir,2)=B(2);
      yfit=XX*B;
      YF(iir,iir:ny)=yfit;
%      dB=abs(B(2)-BB(iir-1));
      if B(2)<=Eps, 
	iR=iir; 
	break; 
      end;
    end
    
    if isempty(iR) % no slope is less than Eps
      iR=find(BB(:,2)==min(BB(:,2)));
    end
    
% Plot time series of annual concentration 
% at the individual location (ik)
% also show linear trends before the cutoff/stauration year
    f_plt=0;
    if f_plt==1
      XY=[1993:2016];
      figure(12); clf;
      axes('Position',[0.09 0.4 0.85 0.52]);
      bar(XY,ava);
      hold
      for iir=1:ny-1
	yf=YF(iir,:);
	yf(yf==0)=nan;
	I=find(~isnan(yf));
	if isempty(I), continue; end
	plot(XY,yf,'r-');
      end
      set(gca,'tickdir','out',...
	      'xtick',[1990:2:2016],...
	      'Fontsize',12);
      stl=sprintf('Pnt %i km, Saturation Year: %i, V.Lev %i',...
		  round(Dst(ik)),iR+xyr(1)-1,ilv);
      title(stl);
      bottom_text(btx,'pwd',1,'Position',[0.01 0.3 0.5 0.05]);
    end
    
    YRTR(ik)=iR+xyr(1)-1; % year when tr = quasi steady state
  end
  
% Fit a polynomial, order=1, 2, 3
  Nyr=YRTR-xyr(1)+1; % # of year to travel distance Dst
  P = polyfit(Nyr,Dst,1); 
  pfit=polyval(P,Nyr);
  PF(ilv,1,:)=pfit;
  
  P = polyfit(Nyr,Dst,2); 
  pfit=polyval(P,Nyr);
  PF(ilv,2,:)=pfit;
  
  P = polyfit(Nyr,Dst,3); 
  pfit=polyval(P,Nyr);
  PF(ilv,3,:)=pfit;
  YearSat(ilv,:)=YRTR; % Saturation years
  
%  figure(20); clf;
%  plot(Nyr,Dst); 
%  hold
%  plot(Nyr,pfit,'r');
%  keyboard
  
  figure(ilv); clf;
  axes('Position',[0.09 0.4 0.85 0.53]);
  pcolor(Dst,xyr,laa); shading flat;
  caxis([-4 1]);
  hold on;
  plot(Dst,YRTR,'r-');
  
  hc=colorbar('southoutside');
  set(hc,'Position',[0.09 0.31 0.85 0.02],...
	 'Ticklength',0.02,'Fontsize',12);
  
  set(gca,'tickdir','out',...
	  'xtick',[0:1000:11000],...
	  'ytick',[1990:2:2018]);
  
  if ilv==1
    lvl='0-50m';
  elseif ilv==2
    lvl='50-150m';
  elseif ilv==3
    lvl='0-bottom';
  end
  
    
  stl=sprintf('log(Ctr), N.Atlantic GSA Route, %s',lvl);
  title(stl);
  bottom_text(btx,'pwd',1,'position',[0.03 0.2 0.6 0.05]);
  
end

% Plot travel speed 
% estimated from the tracer propagation
% fitting a polynomial
% similar figure in
% Belkin 2004, fig. 3
figure(6); clf;
axes('position',[0.2 0.2 0.4 0.72]);
hold on
%CLR=[0 0 0; 0.5 0.5 0.5; 0.8 0.8 0.8];
CLR=[0 0 1; 1 0 0; 0. 0. 0.];

for ilv=1:1
  YRTR=YearSat(ilv,:);
  Nyr=YRTR-xyr(1)+1; % # of year to travel distance Dst
  for ip=1:3
    pf=squeeze(PF(ilv,ip,:));
    clr=CLR(ilv,:);
    plot(Nyr,pf,'Color',clr,'Linewidth',1.6);
    
  end
end
set(gca,'tickdir','out',...
	'xtick',[1:24],...
	'xlim',[1 18],...
	'ylim',[0 11500],...
	'ytick',[0:1000:11000],...
	'xgrid','on',...
	'ygrid','on',...
	'Fontsize',14);
stl=sprintf('Propagation Rate');
title(stl);
bottom_text(btx,'pwd',1,'position',[0.03 0.1 0.6 0.05]);

