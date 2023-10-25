% Statistics of the MLD/ILD: 
% time series of mean MLD 
% over several regions:
% Labrador Sea
% Irmiger Sea
% Greenland 
% Iceland Sea
% BG
% Eurasian Basin
%
% MLD/ILD derived in mld_ArcticOcean.m
% using definition of Mxied Layer
% from Kara et al., 2000
%
% Use monthly mean fields
% Note that original global HYCOM GLBb0.08 outputs 
% have been "remapped" onto ARCc0.08 grid
% see ../REMAP_ARCc/remap_gridGLB2ARC_archv
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

s_fig  = 0; % =0 do not save figure; =1: save figure
yr1 = 1993;
yr2 = 2016;
rg=9806;  % convert pressure to depth, m
dT = 0.5;   % T threshold to calculate d(rho) for MLD - Kara et al., 2003 & 2000
          % dT = 0.2, 0.5
          % 0 - my MLD definition, DD
	  % note for DD method - no ILD calculation, only MLD

pthbin  = '/Net/mars/ddmitry/hycom/GLBb2ARCc0.08/mnth_mean/'; 
pthmat  = '/Net/mars/ddmitry/hycom/GLBb0.08/data_mld/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/'; 
pthfig  = '/Net/mars/ddmitry/hycom/GLBb0.08/fig_MLD/';

TV = '07';
ftopo = sprintf('%sdepth_ARCc0.08_%s.nc',pthtopo,TV); % 
HH  = nc_varget(ftopo,'Bathymetry'); 
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(HH); 
[XX,YY] = meshgrid((1:nn),(1:mm));

% Define indices of the select regions:
% See anls_Greenland sub_define_boxes
BX = sub_define_boxes(HH,LON,LAT,0);

nb = max(size(BX));
for ibx=1:nb
  IJ=BX(ibx).IJ;
  inp = inpolygon(XX,YY,IJ(:,1),IJ(:,2));
  IN = find(inp==1);
  BX(ibx).IN = IN;
end

cc=0;
for year=yr1:yr2
%  fmat = sprintf('%sHYCOManls_monthlyMLD_Arctic_%i',pthmat,year);
  if dT > 0
    fmat = sprintf('%sHYCOManls_monthlyMLD_dT%2.2i_%i',pthmat,dT*10,year);
  else
    fmat = sprintf('%sHYCOManls_monthlyMLD_DD_%i',pthmat,year);
  end    
  fprintf('Loading %s\n',fmat);
  load(fmat);
  
  for im=1:12
    dnmb=MILD(im).TM;
    dv=datevec(dnmb);
    
    MLD = MILD(im).MLD;
    if isfield(MILD,'ILD');
      ILD = MILD(im).ILD;
    else 
      ILD = MLD*0;
    end
    tm  = MILD(im).TM;
    cc=cc+1;
    
    if isempty(MLD),
      fprintf('====  Missing records: %i/%2.2i =======\n',year,im);
      TM(cc,1) = datenum(year,im,1);
      for ibx=1:nb
        mldb(cc,ibx) = nan;
        ildb(cc,ibx) = nan;
      end
      continue
    end
    
    fprintf('Analyzing MLD for %i/%i\n',dv(1),dv(2));

    TM(cc,1)=tm;    
    for ibx=1:nb
      IN=BX(ibx).IN;
      dmm = nanmean(MLD(IN));
      mldb(cc,ibx)=dmm; 
      dmm = nanmean(ILD(IN));
      ildb(cc,ibx)=dmm;
    end
    
    
  end
end

%dT = MILD(1).dT;
CLR=[0,0,0; 0,1,0; 1,0,0; 0,0,1; 1,.6,0; 0 .8 1];
for ibx=1:nb
  NM{ibx}=sprintf('%i %s',ibx,BX(ibx).Name(1:5));
end

DV = datevec(TM);
yr1 = DV(1,1);
yr2 = DV(end,1);
yrs = [0:(yr2-yr1+1)*12-1]/12+yr1;
xlbl=[];
cll=0;
for yr=yrs(1):yrs(end);
  cll=cll+1;
  if mod(cll-1,5)==0
    xlbl{cll}=sprintf('%i',yr);
  else
    xlbl{cll}=' ';
  end
end

% 
% Plot individual regions:
figure(1); clf;
for ibx = 1:nb
  mld = mldb(:,ibx);
  if dT>0
    ttl = sprintf('MLD %s, dT=%3.2f',NM{ibx},dT);
  else
    ttl = sprintf('MLD %s, DukhAlg',NM{ibx});
  end
  
  subplot(3,2,ibx);
  plot(yrs,mld);
  yl1=1.1*min(mld);
  yl2=0;
  title(ttl,'Fontsize',11);
  set(gca,'tickdir','out',...
	  'xlim',[min(yrs)-0.1 max(yrs)+0.1],...
	  'ylim',[yl1 yl2],...
	  'xtick',[yr1:yr2],...
	  'xticklabel',xlbl,...
          'xminortick','on',...
	  'xgrid','on',...
	  'ygrid','on');
end

btx = 'anls_mld_hycom.m';
bottom_text(btx,'pwd',1);
if s_fig>0
  if dT > 0
    fgnm = sprintf('%sGLBb008_monthly_MLDindivid_dt%2.2i',pthfig,dT*10);
  else
    fgnm = sprintf('%sGLBb008_monthly_MLDindivid_DD',pthfig);
  end
  
  fprintf('Saving fig %s\n',fgnm);
  print('-dpng','-r250',fgnm);
end




% Plot monthly climatology
f_allclim=1;
if f_allclim>0
  figure(2); clf;
  axes('position',[0.08 0.55 0.85 0.4]);
  hold on;
  for ibx=1:nb
    clr=CLR(ibx,:);
    plot(mldb(:,ibx),'Linewidth',2,'Color',clr);
  end
  yl1=1.1*min(min(mldb));
  yl2=0;
  set(gca,'xlim',[0.9 12.1],...
	  'ylim',[yl1 0],...
	  'xtick',[1:12],...
	  'tickdir','out');
  if dT>0
    stt=sprintf('HYCOM GLBb0.08, MLD climatology, dT=%3.1f',dT);
  else
    stt=sprintf('HYCOM GLBb0.08, MLD climatology, DD alg.');
  end    
  title(stt,'Fontsize',12);

  axes('position',[0.7 0.65 0.2 0.1]);
  hold on
  for ibx=1:nb
    y1=-ibx;
    clr=CLR(ibx,:);
    plot([0 1],[y1 y1],'-','Linewidth',2,'Color',clr);
    text(2,y1,NM{ibx});
  end
  set(gca,'xlim',[0 6],...
	  'ylim',[-nb-1 0],...
	  'xtick',[],...
	  'ytick',[],...
	  'visible','off');

  if abs(min(min(ildb)))>0
    axes('position',[0.08 0.08 0.85 0.4]);
    hold on;
    for ibx=1:nb
      clr=CLR(ibx,:);
      plot(ildb(:,ibx),'Linewidth',2,'Color',clr);
    end
    yl1=1.1*min(min(ildb));
    yl2=0;
    set(gca,'xlim',[0.9 12.1],...
	    'ylim',[yl1 0],...
	    'xtick',[1:12],...
	    'tickdir','out');

    stt=sprintf('ILD climatology, dT=%3.1f',dT);

    title(stt,'Fontsize',11);
  end
  
  btx='anls_mld_hycom.m';
  bottom_text(btx,'pwd',1);

  if s_fig>0
    fgnm = sprintf('%sGLBb008_monthly_MLDboxes_dt%2.2i',pthfig,dT*10);
    fprintf('Saving fig %s\n',fgnm);
    print('-dpng','-r250',fgnm);
  end
end
