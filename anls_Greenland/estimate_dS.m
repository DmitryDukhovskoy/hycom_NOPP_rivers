% Estimate S change if
% all anomaly is evenly distributed over N.Atl
% in the levels 0-50, 150, 300, 500 m
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers
startup;

close all
clear

s_fig = 0;

nbx = 5; % # of regions
regn = 'ARCc0.08';
expt = 110;  
pthfig  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/fig_NAtl/',expt);
pthsav  = sprintf('/Net/tholia/ddmitry/hycom/ARCc0.08/data_mat%3.3i/',expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
fmout   = sprintf('%sfraction_trcr_NAtl.mat',pthmat);

ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

% Vertical layers
LRS = load('LRS.dat');
nlrs= length(LRS); 

% Grid cell spacing
[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2

% Region:
iGR=[         743        1104
         888        1076
         961         975
        1088         939
	1188         831
        1249         603
        1210         13
         200          20
          70         371
          68        1048
         320         915	      
         465         968
         577        1076];

[XX,YY] = meshgrid((1:nn),(1:mm));
INP = inpolygon(XX,YY,iGR(:,1),iGR(:,2));
IN = find(INP==1 & HH<0);
RG(1).Name = 'SPNA+NordicSeas';
RG(1).IN = IN;
RG(1).iGR= iGR;

% No Nordic Seas
iGR = [         743        1101
         888        1076
         961         975
         992         862
         902         652
         849         582
         911         497
        1086         357
        1249         603
        1210         13
         200          20
          70         371
          68        1048
         320         915	      
         465         968
         577        1076];

INP = inpolygon(XX,YY,iGR(:,1),iGR(:,2));
IN = find(INP==1 & HH<0);
RG(2).Name = 'SPNA no NordicSeas';
RG(2).IN = IN;
RG(2).iGR= iGR;

f_rg=0;
if f_rg==1
  figure(1); clf;
  contour(HH,[0 0],'k');
  hold
  contour(HH,[-1000 -1000],'b');
  iGR=RG(1).iGR;
  plot(iGR(:,1),iGR(:,2),'.-','Linewidth',2);
  iGR=RG(2).iGR;
  plot(iGR(:,1),iGR(:,2),'r.-','Linewidth',1);
  IN = RG(2).IN;
  plot(XX(IN),YY(IN),'c.');
end

  

%[TM,dGR]=sub_read_Greenland_v3; % Greenland anom, km3/yr
frv = sprintf('%sGreenland_cumFWFlux.mat',pthmat);
f_griv = 0; %=1 - rederive cumulative Greenland FWFlux from Bamber
if f_griv==1
  [TMg,Fgr]=sub_read_Greenland_v3; % Greenland total FWF, km3/mo
  % Calculate annual anomalies
  % and cumulative up to date
  DVg = datevec(TMg);
  ngr = length(TMg);
  nyr =ngr/12;
  dmm = reshape(Fgr,[12,nyr]);

  Fyr = sum(dmm); % km3/yr
  Ygr = [DVg(1,1):DVg(end,1)];
  ii=find(Ygr==1992);
  Fmn = mean(Fyr(1:ii));
  ism = find(Ygr==dv0(1,1));
  cFWF = cumsum(Fyr-Fmn);
  fwf0 = cFWF(ism); % km3
  save(frv,'cFWF','Ygr');
else
  fprintf('f_griv %i, Loading %s\n',f_griv,frv);
  load(frv);
%  ii=find(Ygr==1990);
%  ism = find(Ygr==dv0(1,1));
%  fwf0 = cFWF(ism); % km3
end  

ism = find(Ygr==1993);

% Annual mean
fsout = sprintf('%sarc08_%3.3i_mnthS_lrs_%4.4i%2.2i.mat',pthsav,expt,1993,7);
fprintf('Loading S averaged %s\n',fsout);
load(fsout);


ism  = find(Ygr==2016);
fwf0 = cFWF(ism); % km3, note GrFWF anom may be <0 during early 1990s

for irg=1:2
  IN = RG(irg).IN;
  for ilv=1:nlrs-1
% average S over all depths to ilv
    Sav=HH*0;
    for kk=1:ilv
      dZ = abs(LRS(kk,2)-LRS(kk,1));
      Sav = Sav+meanS(kk).Savrg*dZ;
      Sav(Sav==0)=nan;
    end
    hZ  = abs(LRS(ilv,2));
    Sav = Sav./hZ;
    Sav = Sav(IN);
    rr = 1/length(IN); % equally distribute Greenland FWFlx
    Vfw = fwf0*rr*1e9; % km3->m3, 
    
    Vol = Acell(IN).*hZ;
    dS  = -(Sav-(Vol-Vfw).*Sav./Vol); % freshening -> should be negative
    dS(abs(dS)<1e-15)=nan;
    dS(dS>0) = 0; % due to fwf0<0 during early 1990s gives +dS
    
    DltS(irg,ilv)=nanmean(dS); % average dlt S for the regions / depths
  end
end

% Also plot Fram Strait anomaly FW flux 2007-2009
% estimated +1430 of FW surplus flux
% which is ~30% of total Greenland FW flux
% anom estimate for Fram - see plot_FramFW_precip
FrAnom1 = 383; % km3 - estimated FW flux anomaly, in 2007-2008
rfr1 = FrAnom1/fwf0;
FrAnom2= 2511; % km3 anom Fr flux, 2010-2011
rfr2 = FrAnom2/fwf0;
FrAnom3= 3260; % km3 anom Fr Flux, by end 2015
rfr3 = FrAnom3/fwf0;

% Beaufort Gyre anomaly accumulated since 2003 - 2016
BGanom=6600;
rbg = BGanom/fwf0; 

nl=4;
cmp = colormap_gray(nl);
figure(2); clf;
axes('Position',[0.12 0.43 0.7 0.4]);
plot([-2 2],[0 0],'k-');
hold on

dx=0.2;
for irg=1:4
  for ik=1:nl
    switch(irg)
     case(1)
      y0=DltS(1,ik);
      xm=nl/2*dx-irg;
      x1=(ik-1)*dx-xm;
      x2=ik*dx-xm;
      clr=cmp(ik,:);
      aa=patch([x1 x2 x2 x1],[y0 y0 0 0],clr);
     case(2)
% Plot Fram FW anomaly 2007-2008
%    plot([x1 x2],[y0*rfr y0*rfr],'r:');
% Plot Fram FW anomaly 2010-2011
      y0=DltS(1,ik)*rfr2;
      xm=nl/2*dx-irg;
      x1=(ik-1)*dx-xm;
      x2=ik*dx-xm;
      clr=cmp(ik,:);
      aa=patch([x1 x2 x2 x1],[y0 y0 0 0],clr);
     case(3)
% Plot Fram FW anomaly 2015
      y0=DltS(1,ik)*rfr3;
      xm=nl/2*dx-irg;
      x1=(ik-1)*dx-xm;
      x2=ik*dx-xm;
      clr=cmp(ik,:);
      aa=patch([x1 x2 x2 x1],[y0 y0 0 0],clr);
%    plot([x1 x2],[y0*rfr3 y0*rfr3],'m:');
     case(4)
% Plot BG FW anomaly
      y0=DltS(1,ik)*rbg;
      xm=nl/2*dx-irg;
      x1=(ik-1)*dx-xm;
      x2=ik*dx-xm;
      clr=cmp(ik,:);
      aa=patch([x1 x2 x2 x1],[y0 y0 0 0],clr);
%    plot([x1 x2],[y0*rbg y0*rbg],'b:');
    
    end
  end
end


yl1 = 1.1*rbg*min(min(DltS));
xl1 = -(nl/2*dx-1+dx);
xl2 = nl*dx-xm+dx;

sll=sprintf('dS, evenly distrib. GrFW in layers, 2016, Red - Fram dFWFlx, blue - BeaufGyre');
title(sll);
xlb={'GFWA','Fr2010','Fr2015','BG'}; % whole subpolar NA, without Nord.
set(gca,'xticklabel',xlb,...
	'tickdir','out',...
	'ylim',[yl1 0],...
	'xlim',[xl1 xl2],...
	'xtick',[1:4],...
	'ytick',[-0.5:0.05:0],...
	'Fontsize',14);

ptx{1}=sprintf('Fram Anom 2007-2008 %6.1f km3/yr',FrAnom1);
ptx{2}=sprintf('Fram Anom 2010-2011 %6.1f km3/yr',FrAnom2);
ptx{3}=sprintf('Fram Anom 2015 %6.1f km3/yr',FrAnom3);
ptx{4}=sprintf('BG Anom 1993-2016 %6.1f km3/yr',BGanom);

text(0.08, -0.55,ptx,'Fontsize',14);


%lg=legend('0-50','50-150','150-300','300-500');
lvl={'0-50';'50-150';'150-300';'300-500'};
axes('Position',[0.84 0.48 0.125 0.15]);
hold
for ik=1:nl
  x1=0;
  x2=1*dx;
  y1=(ik+0.1)-dx;
  y2=(ik+0.1)+dx;
  clr=cmp(ik,:);
  aa=patch([x1 x2 x2 x1],[y1 y1 y2 y2],clr);
  text(x2+0.25*dx,(y1+y2)/2,lvl(ik),'Fontsize',12);
end
set(gca,'xtick',[],...
	'ytick',[],...
	'xlim',[-0.05 x2+2*dx],...
	'ylim',[1-1.5*dx nl+3*dx],...
	'box','on');

%hc=colorbar('horizontal');
%set(hc,'TickLength',0.018);

btx = 'estimate_dS.m';
bottom_text(btx,'pwd',1,'position',[0.02 0.35 0.4 0.05])

if s_fig==1
  fgnm=sprintf('%sdS_distribFW_SPNA',pthfig);
  fprintf('Saving %s\n',fgnm);
  print('-depsc2',fgnm);
end

