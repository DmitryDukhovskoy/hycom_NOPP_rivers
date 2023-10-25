% Plot Tracer fraction in the N.Atl
% monthly-mean, depth-integrated Tracer mass (by layers)
%  Plotting is done by layers: ilv
%
% For specified locations: same as in Yashayev's time series
%
% Here, use monthly mean depth-integrated within specifed layers LRS
% tracer (mass). Since integration is done exactly 
% within the layers (0-50, 50-150m etc.)
% the concentration is easy to get as dz is fixed
% Ctr = Mass/(dz*Area)
%
% The Tracer Data are extracted in extr_MassTrcr_mnth.m
% mat files: MassTr01_lrs_*.mat for Greenland tracer


addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers
startup;

close all
clear

nTr = 1; 
yr1 = 1993;
yr2 = 2016;
nyrs = yr2-yr1+1;
nbx = 100; % plot boxes =1,..., nbx

s_fig = 0;  % save eps format 
s_mat = 2; % =0 - do not save mat file
           % = 1 save
	   % = 2 skip extraction, only plot output
%dref = datenum(2006,1,15); % reference year/month, = 0 - no reference 
%ilev = 1; %=1 - mixed layer, 0-50; =2 - 50-150m

%nTr = 1; 
% For pfg1 only specify levels:
ilv = 1; % 0-50m
%ilv = 2; % 50-150m
%ilv = 3; % 150-300m
%ilv = 4; % 300-500 m
%ilv = 5; % whole depth <- do not use this
IntSrf = 0; % =1 - integrate from surf to level ilv

% Select figures to plot =0 or =1:
pann = 1; % =0 - plot filtered monthly, =1 - annual mean
pfg0 = 1; % annual time seeer of dS fooor 0-50,50-200, 50-500 sim to Yash
pfg1 = 0; % bar diagrams FWC per dZrf m of water column
pfg2 = 0; % time ser of FWC cont change for all levels (mat files have to be created)
pfg3 = 0; % time ser. dlt S change within the layer, all layers


% Vertical layers
LRS = load('LRS.dat');
nlrs= length(LRS); 


zz1 = LRS(ilv,1);
zz2 = LRS(ilv,2);

dz=abs(zz2-zz1);

fprintf('Time Ser. Tracer Mass budget, ilv=%i, %i - %i, nTr=%i\n',ilv,zz1,zz2,nTr);


regn = 'ARCc0.08';
expt = 110;  
pthfig  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/fig_NAtl/',expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
pthsav  = sprintf('/Net/tholia/ddmitry/hycom/ARCc0.08/data_mat%3.3i/',expt);
fll     = sprintf('%sGrTrFrct_GrWVol_YashRegions_Sref1993_',pthmat);
%fll   = sprintf('%sGrTrFrct_GrWVol_YashRegions_NAtl_',pthmat);
fmout   = sprintf('%slev%2.2i.mat',fll,ilv);
btx     = 'fraction_Yash_MassTrcr_NAtl.m';

LRS = load('LRS.dat');
nlrs= length(LRS); 

dZrf = 50; 
if s_mat<2

  ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
  HH  = nc_varget(ftopo,'Bathymetry');
  LON = nc_varget(ftopo,'Longitude');
  LAT = nc_varget(ftopo,'Latitude');
  [mm,nn]=size(LON);

  % Grid cell spacing
  [DX,DY]=sub_dx_dy(LON,LAT);
  Acell=DX.*DY; % Grid cell area, m2


  xlim1 = 20;
  xlim2 = nn;
  ylim1 = 5;
  ylim2 = 2000;

  hmsk=HH;
  hmsk(HH<0)=nan;

  % Define Regions:
  f_pltbox=0;
  BX = sub_define_Yash_reg(HH,LON,LAT,f_pltbox);
  if isempty(nbx) | nbx>20, nbx=length(BX); end

  [XX,YY] = meshgrid((1:nn),(1:mm));
  for ib=1:nbx
    iBG = BX(ib).IJ;
    INP = inpolygon(XX,YY,iBG(:,1),iBG(:,2));
    IN = find(INP==1);
    BX(ib).IN = IN;
  end


  chck=0;
  if chck>0
    figure(10); clf;
    contour(HH,[0 0],'k');
    hold on;
    contour(HH,[-5000:1000:-10],'Color',[0.5 0.5 0.5]);
  %  contour(HH,[-1000 -1000],'b');
    for ib=1:nbx
      iBG=BX(ib).IJ;
      plot(iBG(:,1),iBG(:,2),'r.-');
    end
  %  contour(LON,[-153 -153],'g');
    contour(LON,[-150:30:150],'Color',[0.8 0.8 0.8]);
  %  contour(LAT,[70 70],'g');
    contour(LAT,[50:10:80],'Color',[0.9 0.9 0.9]);
    [JJ,II]=ind2sub(size(HH),IN);

    for ib=1:nbx
      X=BX(ib).XY(:,1);
      Y=BX(ib).XY(:,2);
      iBG=BX(ib).IJ;
      x0=mean(X);
      y0=mean(Y);
      for kk=1:4
	x=X(kk);
	y=Y(kk);
	dx=x-x0;
	dy=y-y0;

	i0=iBG(kk,1)+20*sign(dx);
	j0=iBG(kk,2)+20*sign(dy);
	stl{1}=sprintf('%5.2fN',y);
	stl{2}=sprintf('%5.2fE',x);
	tx=text(i0,j0,stl,'Fontsize',12);
	set(tx,'Color',[0 0 1]);
      end

    end

    axis('equal');

  end

  % Read in cumulative Greenland FW flux anomaly
  % mat file created in dSregions_map_v2.m
  frv = sprintf('%sGreenland_cumFWFlux.mat',pthmat);
  fprintf('Loading %s\n',frv);
  load(frv);

% ----------------------------------
% Create time series of Tracer Mass
% and Greenland FW present in the
% specified regions
% ----------------------------------
% 
  cc=0;
%  clear sumT
  for iyr=yr1:yr2
    for imo=1:12
%pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
% fmat = sprintf('%sMassTr%2.2i_lrs_%i%2.2i.mat',pthmat,nTr,iyr,imo);
      dnmb = datenum(iyr,imo,15);
      dv0  = datevec(dnmb);
      yr   = dv0(1);
      iday = dnmb-datenum(yr,1,1)+1;
      rr = sub_fraction_tracerMass(expt,regn,Acell,HH,dnmb,nTr,ilv,dz,IntSrf);
      rr(isnan(rr))=0;
      
      if isempty(rr)
	cc=cc+1;
        for ib=1:nbx % regions
	  BX(ib).GrFWcont_m_dZrf(cc,1) = nan;
	end
	continue;
      end

%
% Get layer-averaged S - as a reference S
% use 1 year (1993) 
% to calculate dlt S
% prepared in anls_TS/mnthly_arc08_layers_S.m
      YR=dv0(1);
      mo=dv0(2);
%      fsout = sprintf('%sarc08_%3.3i_mnthS_lrs_%4.4i%2.2i.mat',pthsav,expt,YR,mo);
      fsout = sprintf('%sarc08_%3.3i_mnthS_lrs_%4.4i%2.2i.mat',pthsav,expt,1993,mo);
      fprintf('Loading S averaged %s\n',fsout);
      load(fsout);
      Sav = meanS(ilv).Savrg;
      Sav(Sav==0)=nan;
      hZ = abs(LRS(ilv,2)-LRS(ilv,1));

% Greenland runoff anomaly for given date
%  ii=find(Ygr==1990);
      ism = find(Ygr==dv0(1,1));
% interpolate into months:
      df=(cFWF(ism)-cFWF(ism-1))/12; 
      fwf0 = cFWF(ism-1)+df*imo; % km3, note GrFWF anom may be <0 during early 1990s
      

% Estimate volume of Greenland surplus FW in grid cell
% = ratio of FW tracer*total FWF (integrated over time 1990-date requested)
      Vfw = fwf0*rr*1e9; % m3 

% Estiamte dlt(S) due to Greenland FWFlux anomaly
      Vol = Acell.*hZ;
      dS=-(Sav-(Vol-Vfw).*Sav./Vol); % freshening -> should be negative
      dS(abs(dS)<1e-15)=nan;
      dS(dS>0)=0; % due to fwf0<0 during early 1990s gives +dS
      
      
      cc=cc+1;
      for ib=1:nbx % regions
	IN=BX(ib).IN;
	Rarea = sum(Acell(IN));   % total area region
	Fgr   = sum(rr(IN)*fwf0*1e9); % m3, vol of Gr. FW in the region
	FgrM  = (Fgr/Rarea)/dz*dZrf; % m of Gr. FW in dzRf m of water
	BX(ib).dZrf = dZrf;  % reference layer thickness for FWC normalization
	BX(ib).GrFWcont_m_dZrf(cc,1) = FgrM; % avrg Tr content in region
	BX(ib).dS(cc,1) = nanmean(dS(IN)); % mean dlt S in the region
	BX(ib).dSm(cc,1)= min(dS(IN)); % min dlt S in the region
%	if nanmean(dS(IN))>1e-4, fprintf('dS>0\n'), keyboard; end
%	BX(ib).TM(cc,1)=dnmb;
      end

    end % month
  end   % year
  if s_mat==1
    fprintf('Saving %s\n',fmout);
    save(fmout,'BX');
  end
end

if s_mat==1; return; end;

clear ABX
for ilv=1:4
%  fmout   = sprintf('%sGrTrFrct_GrWVol_YashRegions_NAtl_lev%2.2i.mat',pthmat,ilv);
  fmout   = sprintf('%slev%2.2i.mat',fll,ilv);
  fprintf('Loading %s\n',fmout);
  load(fmout);
  if isempty(nbx) | nbx>20, nbx=length(BX); end
  for ibx=1:nbx
    nm = BX(ibx).Name;
    ABX(ibx).Name = nm;
    grw = BX(ibx).GrFWcont_m_dZrf;
    ABX(ibx).GrFWCm(ilv,:)=grw;
    dS = BX(ibx).dS;
    ABX(ibx).dS(ilv,:) = dS;
  end
end

% Prepare time series 0-50, 50-200, 
% 50-500 m
for ibx=1:nbx
  ilv=5;
  dS = ABX(ibx).dS;
  ds1=100/150*dS(2,:)+50/150*(dS(3,:)); % dS in 50-200 m
  ABX(ibx).dS(ilv,:)=ds1;
  
  ilv=6;
  ds1=100/450*dS(2,:)+150/450*dS(3,:)+200/250*dS(4,:); % dS in 50-500 m
  ABX(ibx).dS(ilv,:)=ds1;
end      

CLR1 = [0 0.4 0.6; ...
       0.7 0.3 0; ...
       0. 0.7 0.5];
  
SITES=struct; % for Igor Yash. 
nst=0;
pyrs=[yr1:yr2]+0.5;
if pfg0==1   % plot 0-50, 50-200 m,   50-500m
  for ib=1:nbx
    if ib==12, continue; end; % skip repition of SW Iceland Shelf
    nm = BX(ib).Name;
    figure(ib); clf;
    axes('Position',[0.09 0.45 0.85 0.45]);
    hold on;
    yl1=0;
    cl=0;
    
    nst=nst+1;
    SITES(nst).Name=nm;
    cLRS={'0-50','50-200','50-500'};
    
    for ilv=[1,5,6]
      cl=cl+1;
      dS = ABX(ib).dS(ilv,:);
      inn = find(~isnan(dS));
      if ~isempty(inn),
        dS(isnan(dS))=0;
      end
      
      clr = CLR1(cl,:);
      A = reshape(dS,[12,nyrs]);
      dmm = mean(A);
      yrs = [yr1:yr2]+0.5;
      dmm(dmm>0)=0;
      
% Fit linear polynom over segment
      yy=dmm';
      xx=[1:length(yy)];
      ons=xx*0+1;
      X=[ons',xx'];
      B = regress(yy,X);
      yft = X*B;
      xxy = xx'+yrs(1)-1;
      
      pp=plot(yrs,dmm,'Linewidth',2,'Color',clr);
      plg(cl)=pp;
      plot(yrs,dmm,'k.','Color',clr,'Markersize',18);
%      clrL=clr+[.2 .2 .2];
      clrL=clr;
      plot(xxy,yft,'Color',clrL,'Linewidth',3);
      mdd=min(dmm);
      yl1=min([yl1,mdd]);
      txt{cl}=sprintf('%i, dS 10yr %6.5d',cl,B(2)*10);

      SITES(nst).Layers(cl)=cLRS(cl);
      SITES(nst).dS_monthly(cl,:)=dS;
      SITES(nst).dS_annual(cl,:)=dmm;
      SITES(nst).Years=yrs;
      SITES(nst).Regress_cfnt(cl,:)=B';
      
    end

    set(gca,'tickdir','out',...
	    'xlim',[1993 2017],...
	    'xtick',[1993:2017],...
	    'ylim',[1.05*yl1 0],...
	    'Fontsize',15);
    
    if mod(ib,2)==0,
      set(gca,'YAxisLocation','right');
    end
    

    hll=legend(plg,'0-50','50-200','50-500',...
	   'Location','SouthOutside');
    set(hll,'Position',[0.70 0.24 0.1 0.1],'Fontsize',12);

    stl=sprintf('%s, dS Est. from Tracer Budget',nm);
    title(stl);
    
    axes('Position',[0.05 0.15 0.2 0.1]);
    text(1,1,txt);
    set(gca,'xlim',[1 1.25],'visible','off','Fontsize',14);
    set(gcf,'Position',[654 637 1655 645]);

    bottom_text(btx,'pwd',1,'position',[0.02 0.3 0.4 0.1]);
%keyboard

    if s_fig==1
      fgnm=sprintf('%sdS_GrExpt_YashPnts_rg%2.2i',pthfig,ib);
      fprintf('Saving fig %s\n',fgnm);
      print('-depsc2',fgnm);
    end
    

  end  % regions

  % Save time series for I. Yashaev
%    fyash=sprintf('%shycom_dS_regions.mat',pthmat);
%    fprintf('Saving %s\n',fyash);
%    save(fyash,'SITES');
    

end

tmm=[yr1:yr2];
if pfg1==1
  for ibx=1:nbx
    nm = ABX(ibx).Name;
    figure(ibx); clf;
    axes('Position',[0.08 0.55 0.88 0.35]);
    grw = ABX(ibx).GrFWCm(ilv,:);
    I=find(grw<1e-20);
    grw(I)=nan;
    A=reshape(grw,12,24);
    A=nanmean(A,1);
    bar(tmm,A,0.97,'Facecolor',[0.7 0.7 0.7],'edgecolor','none');
    
    set(gca,'tickdir','out',...
	    'xlim',[1992.5 2016.5],...
	    'xtick',[1993:2016],...
	    'ylim',[0 1.1*max(A)],...
	    'FontSize',14);
    stl = sprintf('%s, GreenlFWC m per %3.1i m, Lr: %i-%i m',...
		  nm,dZrf,abs(zz1),abs(zz2));
    title(stl);
    
    bottom_text(btx,'pwd',1,'Position',[0.02 0.4 0.5 0.1]);
  end
end


% Butterworth filter
% output freq. 1/12 mo^-1
% 12 mo cutoff = 12/12=1 -> in Matlab Wn=1/6
% 
Wn = 1/6; % cutoff freq 6 mo: 2/6, 1yr=1/6
[Bf,Af] = butter(9,Wn,'low');

CLR = [0 0.4 0.6; ...
       0.7 0.3 0; ...
       0.8 0.7 0.; ...
       0. 0.8 0.4];

% Time ser of FWC change due to Greenl FW Flux anomaly
yrs=[1993:1/12:2016.99];
if pfg2==1
  for ib=1:nbx
    nm = ABX(ib).Name;
    figure(ib+10); clf;
    axes('Position',[0.09 0.45 0.85 0.45]);
    hold on;
    lmx=0;
    for ilv=1:4
      dFWC = ABX(ib).GrFWCm(ilv,:);
      dFWC(dFWC<0)=0;
      dmm = filtfilt(Bf,Af,dFWC);
      lmx=max([lmx, max(dmm)]);
      clr = CLR(ilv,:);
      plot(yrs,dmm,'Linewidth',1.6,'Color',clr);
    end

    set(gca,'tickdir','out',...
	    'xlim',[1993 2017],...
	    'xtick',[1993:2017],...
	    'ylim',[0 1.1*lmx],...
	    'xgrid','on',...
	    'ygrid','on',...
	    'Fontsize',12);

    hll=legend('0-50','50-150','150-300','300-500',...
	   'Location','SouthOutside');
    set(hll,'Position',[0.52 0.26 0.1 0.1],'Fontsize',12);

    stl=sprintf('%s, dFWC from Tracer, per %4.1fm Water',nm,dZrf);
    title(stl);

    bottom_text(btx,'pwd',1,'position',[0.02 0.3 0.4 0.1]);

  end

  
end



% Plot time series of dltS 
% similar to dltS from Greenland on/off comparison
% see anls_TS/dltS_tser_LrAbrg_GreenldExp_NAtlRegions.m
nyrs=24;
if pfg3==1
  for ib=1:nbx
    nm = BX(ib).Name;
    figure(ib+20); clf;
    axes('Position',[0.09 0.45 0.85 0.45]);
    hold on;
    for ilv=1:4
      dS = ABX(ib).dS(ilv,:);
      inn = find(~isnan(dS));
      if ~isempty(inn),
        dS(isnan(dS))=0;
      end
      
      clr = CLR(ilv,:);
      if pann == 0;
        dmm = filtfilt(Bf,Af,dS);
      else
	A = reshape(dS,[12,nyrs]);
	dmm = mean(A);
	yrs = [yr1:yr2];
      end
      dmm(dmm>0)=0;
      pp=plot(yrs,dmm,'Linewidth',1.8,'Color',clr);
      plg(ilv)=pp;
      if pann==1
	plot(yrs,dmm,'k.','Color',clr,'Markersize',15);
      end
      
    end

    set(gca,'tickdir','out',...
	    'xlim',[1993 2017],...
	    'xtick',[1993:2017],...
	    'xgrid','on',...
	    'ygrid','on',...
	    'Fontsize',15);

    hll=legend(plg,'0-50','50-150','150-300','300-500',...
	   'Location','SouthOutside');
    set(hll,'Position',[0.70 0.24 0.1 0.1],'Fontsize',12);

    stl=sprintf('%s, dS Est. from Tracer Budget',nm);
    title(stl);
    set(gcf,'Position',[654 637 1655 645]);

    bottom_text(btx,'pwd',1,'position',[0.02 0.3 0.4 0.1]);


  end
end  