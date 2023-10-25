% Plot Tracer fraction in the N.Atl
% monthly-mean, depth-integrated Tracer mass (by layers)
%  Plotting is done by layers: ilv
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
nbx = 5; % plot boxes =1,..., nbx

s_fig = 0;
s_mat = 2; % =0 - do not save mat file
           % = 1 save
	   % = 2 skip extraction, only plot output
%dref = datenum(2006,1,15); % reference year/month, = 0 - no reference 
%ilev = 1; %=1 - mixed layer, 0-50; =2 - 50-150m

nTr = 1; 
% For pfg1 only specify levels:
ilv = 1; % 0-50m
%ilv = 2; % 50-150m
%ilv = 3; % 150-300m
%ilv = 4; % 300-500 m
%ilv = 5; % whole depth <- do not use this

% Select figures to plot =0 or =1:
pfg1 = 0; % bar diagrams FWC per dZrf m of water column
pfg2 = 0; % time ser of FWC cont change for all levels (mat files have to be created)
pfg3 = 1; % time ser. dlt S change within the layer, all layers


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
fmout   = sprintf('%sGrTrFrct_GrWVol_Regions_NAtl_lev%2.2i.mat',pthmat,ilv);
btx     = 'fraction_MassTrcr_NAtl.m';

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
  BX = sub_define_boxes(HH,LON,LAT,f_pltbox);

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
      rr = sub_fraction_tracerMass(expt,regn,Acell,HH,dnmb,nTr,ilv,dz);
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
% to calculate dlt S
% prepared in anls_TS/mnthly_arc08_layers_S.m
      YR=dv0(1);
      mo=dv0(2);
      fsout = sprintf('%sarc08_%3.3i_mnthS_lrs_%4.4i%2.2i.mat',pthsav,expt,YR,mo);
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
      fwf0 = cFWF(ism-1)+df*imo; % km3
      

% Estimate volume of Greenland surplus FW in grid cell
% = ratio of FW tracer*total FWF (integrated over time 1990-date requested)
      Vfw = fwf0*rr*1e9; % m3 

% Estiamte dlt(S) due to Greenland FWFlux anomaly
      Vol = Acell.*hZ;
      dS=-(Sav-(Vol-Vfw).*Sav./Vol); % freshening -> should be negative
      dS(abs(dS)<1e-20)=nan;
      
      
      cc=cc+1;
      for ib=1:nbx % regions
	IN=BX(ib).IN;
	Rarea = sum(Acell(IN));   % total area region
	Fgr   = sum(rr(IN)*fwf0*1e9); % m3, vol of Gr. FW in the region
	FgrM  = (Fgr/Rarea)/dz*dZrf; % m of Gr. FW in dzRf m of water
	BX(ib).dZrf = dZrf;  % reference layer thickness for FWC normalization
	BX(ib).GrFWcont_m_dZrf(cc,1) = FgrM; % avrg Tr content in region
	BX(ib).dS(cc,1) = nanmean(dS(IN)); % mean dlt S in the region
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
  fmout   = sprintf('%sGrTrFrct_GrWVol_Regions_NAtl_lev%2.2i.mat',pthmat,ilv);
  fprintf('Loading %s\n',fmout);
  load(fmout);
  for ibx=1:nbx
    nm = BX(ibx).Name;
    ABX(ibx).Name = nm;
    grw = BX(ibx).GrFWcont_m_dZrf;
    ABX(ibx).GrFWCm(ilv,:)=grw;
    dS = BX(ibx).dS;
    ABX(ibx).dS(ilv,:) = dS;
  end
end

      

  
tmm = [yr1:yr2];
nyr=24;
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
  for ib=1:5
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

if pfg3==1
  for ib=1:5
    nm = BX(ib).Name;
    figure(ib+20); clf;
    axes('Position',[0.09 0.45 0.85 0.45]);
    hold on;
    for ilv=1:4
      dS = ABX(ib).dS(ilv,:);
      dS(isnan(dS))=0;
      dmm = filtfilt(Bf,Af,dS);
      clr = CLR(ilv,:);
      plot(yrs,dmm,'Linewidth',1.6,'Color',clr);
    end

    set(gca,'tickdir','out',...
	    'xlim',[1993 2017],...
	    'xtick',[1993:2017],...
	    'xgrid','on',...
	    'ygrid','on',...
	    'Fontsize',12);

    hll=legend('0-50','50-150','150-300','300-500',...
	   'Location','SouthOutside');
    set(hll,'Position',[0.52 0.26 0.1 0.1],'Fontsize',12);

    stl=sprintf('%s, dS Est. from Tracer Budget',nm);
    title(stl);

    bottom_text(btx,'pwd',1,'position',[0.02 0.3 0.4 0.1]);

  end
end  