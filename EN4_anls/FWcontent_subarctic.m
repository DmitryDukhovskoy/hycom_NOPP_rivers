% FW content in subarctic seas
% http://hadobs.metoffice.com/en4/index.html
%The EN4 dataset consists of two products:
%Observed subsurface ocean temperature and 
% salinity profiles with data quality information, and,
%Objective analyses formed from the profile data with uncertainty estimates.
%Data are available from 1900 to the present and 
% there are separate files for each month.
% Please read 'Good, S. A., M. J. Martin and N. A. Rayner, 2013. EN4: 
% quality controlled ocean temperature and salinity profiles and 
% monthly objective analyses with uncertainty estimates, 
% Journal of Geophysical Research: Oceans, 118, 6704-6716, 
% doi:10.1002/2013JC009067' for details of how the dataset was constructed.
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

s_fig = 0; 
s_mat = 1; % =0 - load existing mat file
           % =1 calculate FWC
f_month=0; % =1 - plot monthly FWC + annual on top


%pthdat = '/Net/yucatan/tachanat/ocean_analysis/EN4/EN4_extract/';
pthdat = '/Net/kronos/ddmitry/EN4/';
%pthmat = '/Net/ocean/ddmitry/vector_winds/dataM/'; 
pthmat = '/nexsan/people/ddmitry/data_mat/';
pthmat2= '/Net/mars/ddmitry/hycom/GLBb0.08/data_mat/';
pthfig = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/fig_EN4/';
en4v   = 'EN.4.2.0.f.analysis.g10'; % EN4 version
txtb   = 'FWcontent_subarctic.m';

ftopo = sprintf('%sGLBb_T07_subset_Natl.mat',pthmat2);
fprintf('Loading topo %s\n',ftopo);
TH = load(ftopo);
HH = TH.HH;
LONH = TH.LON;
LATH = TH.LAT;

%Sref = 35.5;
%Sref = 35; 
Sref = 34.8; 
yr1=1990;
yr2=2016;
fmat = sprintf('%sFWC_subarctic_EN4_Sref%3.3i.mat',pthmat,round(Sref*10));

fprintf('        Sref = %6.2f\n',Sref);

fnm = sprintf('%s%s.201302.nc',pthdat,en4v);

S = double(nc_varget(fnm,'salinity'));
S = squeeze(S(1,1,:,:));

% Find sections:
LON = nc_varget(fnm,'lon');
LAT = nc_varget(fnm,'lat');
ZM  = -1*nc_varget(fnm,'depth');
lz = length(ZM);
dzm = abs(diff(ZM));
dzm(lz) = dzm(lz-1);

clear ZZ
ZZ(1,1) = 0;
for kk=1:lz
  ZZ(kk+1) = -(abs(ZM(kk))+abs(0.5*dzm(kk)));
end
ZZ=ZZ(:);
DZ=abs(diff(ZZ)); % layer thicknesses

I = find(LON>180);
LON(I)=LON(I)-360;

FWC = sub_regions;
nR = length(FWC);

%keyboard

% Find indices:
[LN,LT] = meshgrid(LON,LAT);
[DX,DY] = sub_dx_dy(LN,LT);
Acell = DX.*DY*1e-6; % km2

for i=1:nR
  XY = FWC(i).XY;
%  if i ~= 4
  pl = inpolygon(LN,LT,XY(:,1),XY(:,2));
  IN = find(pl==1);
  
  if i==4,
    nL=length(LON);
    nT=length(LAT);
    x1=min(XY(:,1));
    x2=max(XY(:,1));
    y1=min(XY(:,2));
    y2=max(XY(:,2));
    i1=max(find(LON<=x1))
    i2=max(find(LON<=x2 & LON>0));
    j1=max(find(LAT<=y1));
    j2=max(find(LAT<=y2));
    I1=[i1:nL];
    I2=[1:i2];
    J1=[j1:j2];
    II=[I1,I2];
    cx=0;
    for ii=II
      for jj=J1
	cx=cx+1;
	I=sub2ind(size(LN),jj,ii);
	IN(cx,1)=I;
      end
    end
  end
  
  FWC(i).IN = IN;
%  else

  FWC(i).ZZ   = ZZ;
  FWC(i).Code = '/hycom_NOPP/EN4_anls/FWCcontent_subarctic.m';
  
% Find HYCOM indices corresponding to selected regions:
  for jj=1:4
    x=FWC(i).XY(jj,1);
    y=FWC(i).XY(jj,2);
    D = distance_spheric_coord(LATH,LONH,y,x);
    [j0,i0] = find(D==min(min(D)));
    FWC(i).IJs(jj,1)=i0;
    FWC(i).IJs(jj,2)=j0;
  end
  FWC(i).IJs(jj+1,1)=FWC(i).IJs(1,1);
  FWC(i).IJs(jj+1,2)=FWC(i).IJs(1,2);
end
%keyboard
% =======================
f_chck=0;
if f_chck>0
  figure(1); clf;
  contour(LONH,LATH,HH,[0 0],'k');
  hold on;
  for i=1:nR
    XY = FWC(i).XY;
    plot(XY(:,1),XY(:,2),'r.-');
    IN = FWC(i).IN;
    plot(LN(IN),LT(IN),'g.');
  end
end

f_reg=0;
if f_reg>0
% Plot regions
  nf=20;
  nrplt=[];
  sub_plot_regions(nf,FWC,HH,LONH,LATH,nrplt);
  bottom_text(txtb,'pwd',1);
  

end

% =======================



cc=0;
for ii=yr1:yr2
  for im=1:12
%    if im>2 & im<11, continue; end
    cc=cc+1;
    YRPLT(cc,1)=ii;
    YRPLT(cc,2)=im;
  end
end
nrc=cc;

if s_mat>0
  cc=0;
  for ik=1:nrc
    YR=YRPLT(ik,1);
    IM=YRPLT(ik,2);
    dnmb = datenum(YR,IM,1);
    fnm = sprintf('%s%s.%4.4i%2.2i.nc',pthdat,en4v,YR,IM);

    fprintf('EN4: Reading %i/%i, %s\n',YR,IM,fnm);

    SAL = double(nc_varget(fnm,'salinity'));
    TT = double(nc_varget(fnm,'temperature'))-273.15; % K -> C
    cc=cc+1;

    SAL = squeeze(SAL(1,:,:,:));
    TT = squeeze(TT(1,:,:,:));
    [ll,mm,nn]=size(SAL);

    
    for i=1:nR
      IN = FWC(i).IN;
      S = SAL(:,IN);
      T = TT(:,IN);
      nm = FWC(i).Name;

% Follow ~Haine et al., 2015 definition of FWC
% The better way would be to interpolate
% between S values to find the exact depth
      Fwc = 0; 
      Isrf = [];
      for k = 1:ll-1
	dz = DZ(k);
	ss1=squeeze(S(k,:));
	ss2=squeeze(S(k+1,:));
	ac = Acell(IN);
	Ib1=find(ss1>Sref);
%	if k==1, Isrf=Ib1; end; % exclude region with surf S>Sref
	ss1(Ib1)=nan;
	ss2(Ib1)=nan;
	Ib2=find(ss2>Sref); % CHeck next layer see if S > Sref there
	fwc=dz.*(Sref-ss1)/Sref; % m
  % Interpolate to find depth of Sref
  % Take care about layers where S becomes S>Sref
  % add FWC below layer interface to the depth of Sref
	if ~isempty(Ib2)
	  zz1 = ZZ(k);
	  zz2 = ZZ(k+1);
	  dZz = abs(zz2-zz1);
	  dS = ss2-ss1;
	  zzSref = zz1+dZz./dS.*(Sref-ss1);
	  dltZs  = abs(zz2-zzSref);
	  Smn=0.5*(ss1+Sref); % for integration from 
	  dltFwc = dltZs.*(Sref-Smn)/Sref;
	  fwc(Ib2) = fwc(Ib2)+dltFwc(Ib2);
	end
	fwc(Ib1)=0;    
	fwc=fwc(:);
	fwc(Isrf) = 0;
	Fwc=Fwc+nansum(fwc.*ac)*1e-3; %FWC volume-integrated, km3
%if i==6, keyboard; end;
      end % k level
%      Fwc(Isrf)=0; 
  %  Fwc(Ilnd)=nan; % exclude not needed regions
      fprintf('====> %s: %s, FWC=%8.2f m\n',datestr(dnmb),nm,Fwc/sum(ac)*1e3); 
%keyboard
% Spatial average T/S profile:
      Sav = nanmean(S,2);
      Tav = nanmean(T,2);

      FWC(i).TM(cc,1)      = dnmb;
      FWC(i).Fwc_km3(cc,1) = Fwc;
      FWC(i).Sref          = Sref;
      FWC(i).Area_km2      = sum(ac);
      FWC(i).Fwc_m(cc,1)   = Fwc/sum(ac)*1e3;
      FWC(i).S_spatial_av(:,cc) = Sav;
      FWC(i).T_spatial_av(:,cc) = Tav;
    end; % Regions
%keyboard    
  end % Time

  if s_mat==1
    fprintf('Saving %s\n',fmat);
    save(fmat,'FWC');
  end
  
else
  fprintf('Loading %s\n',fmat);
  load(fmat);
end

TM = FWC(1).TM;
DV = datevec(TM);
YR1 = DV(1,1);
YR2 = DV(end,1);
yrs = [0:length(TM)-1]./12+YR1;


% Calculate yearly mean FWC:
clear Fm
for i=1:nR
  Fwc = FWC(i).Fwc_m;
  a2=length(Fwc);
  Nyr=a2/12;

  b=reshape(Fwc,[12,Nyr]);
  Fm(i,:) = nanmean(b,1);
end



xlbl=[];
cll=0;
for yr=yrs(1):yrs(end);
  cll=cll+1;
  if mod(yr,5)==0
    xlbl{cll}=sprintf('%i',yr);
  else
    xlbl{cll}=' ';
  end
end


if f_month>0
  for ik=1:nR
    figure(ik); clf;
    axes('position',[0.1 0.5 0.85 0.4]);
    ss1=FWC(ik).Fwc_m;
    plot(yrs,ss1,'-');
    hold on;
    for jr=1:Nyr
      x1=YR1+jr-1;
      x2=x1+0.999;
      plot([x1 x2],[Fm(ik,jr) Fm(ik,jr)],'r-');
    end

    set(gca,'xlim',[min(yrs)-0.1 max(yrs)+0.1],...
	    'ylim',[0.999*min(ss1) 1.001*max(ss1)],...
	    'tickdir','out',...
	    'xtick',[yr1:yr2],...
	    'xticklabel',xlbl,...
	    'xminortick','on',...
	    'xgrid','on',...
	    'ygrid','on');

    nm=FWC(ik).Name;
    stl = sprintf('EN4, FWC, %i-%i, %s, Sref=%4.1f',...
		  YR1,YR2,nm,Sref);

    title(stl);
    btx = 'TS_convection_regions.m';
    bottom_text(btx,'pwd',1,'position',[0.08 0.3 0.7 0.2]);

    if s_fig>0
      fgnm=sprintf('%sFWCmonth_EN4_%s',pthfig,nm);
      fprintf('Saving %s\n',fgnm);
      print('-dpng','-r250',fgnm);
    end

  end
end
%
% Plot FWC annual bar plot:
%
YR = [yr1:yr2];
for ik=1:nR
  figure(ik+10); clf;
  axes('position',[0.1 0.5 0.85 0.4]);
  nm=FWC(ik).Name;
  fwv = Fm(ik,:);

  hb=bar(YR,fwv,0.9);
  set(hb,'Facecolor',[0.2 0.5 1]);
  stl = sprintf('EN4, FWC, %s, Sref=%4.1f',...
		  nm,Sref);
  title(stl,'Fontsize',11);
  yl1=floor(min(fwv)-0.1*min(fwv));
  yl2=ceil(max(fwv)+0.1*max(fwv));
  dyl=1;
  if (yl2-yl1)>10,
    dyl=2;
  end
  
%  set(gca,'xlim',[YR(1)-0.5 YR(end)+0.5],...
  set(gca,'xlim',[1992.5 YR(end)+0.5],...
	  'xtick',[YR(1):YR(end)],...
	  'xticklabel',xlbl,...
	  'ylim',[yl1 yl2],...
	  'ytick',[0:dyl:50],...
	  'tickdir','out',...
	  'fontsize',14);
%  ylabel('FWC, km3');

  bottom_text(txtb,'pwd',1,'position',[0.08 0.3 0.7 0.2]);
  
  if s_fig>0
    fgnm=sprintf('%sFWCann_EN4_%s',pthfig,nm);
    fprintf('Saving %s\n',fgnm);
    print('-dpng','-r250',fgnm);
  end

end

sclr = [0,0,0.8; 0,0.8,1; 1,1,1; 1,0.8,0; 0.8,0,0];
c1 = -1;
c2 = 1;
nint = 200;
CMP = create_colormap(nint,c1,c2,sclr);
cmpT = CMP.colormap;
cntT = CMP.intervals;

cs1 = -0.1;
cs2 = 0.1;
nint = 200;
CMP = create_colormap(nint,c1,c2,sclr);
cmpS = CMP.colormap;
cntS = CMP.intervals;


% Plot time vs T/S vertical profiles
DV  = datevec(TM);
TMx = DV(:,1)+DV(:,2)/12;
nrc = length(TMx);
for ik=1:nR
  Tav = FWC(ik).T_spatial_av;
  Sav = FWC(ik).S_spatial_av;
  nm  = FWC(ik).Name;
  Tav = [Tav(1,:);Tav];
  Sav = [Sav(1,:);Sav];
  
% Estimate climatology
  clear Tcl Scl dltT dltS
  for k=1:lz+1
    dmm = reshape(Tav(k,:,:),[12,27]);
    dmm = nanmean(dmm,2);
    Tcl(k,:) = dmm;
    dmm = reshape(Sav(k,:,:),[12,27]);
    dmm = nanmean(dmm,2);
    Scl(k,:) = dmm;
  end
  
  for jc=1:nrc
    mo=mod(jc,12);
    if mo==0, mo=12; end;
    dltT(:,jc) = Tav(:,jc)-Tcl(:,mo);
    dltS(:,jc) = Sav(:,jc)-Scl(:,mo);
  end
  
  figure(ik+1); clf;
% Plot dlt T
  axes('Position',[0.08 0.55 0.85 0.4]);
  pcolor(TMx,ZZ,dltT); shading flat
  caxis([c1 c2]);
  colormap(cmpT);
  CB = colorbar;
  set(CB,'Position',[0.935 0.55 0.01 0.4],...
	 'TickLength',0.04);
  
  txt = sprintf('%s, dltT (base: %i-%i)',nm,floor(TMx(1)),floor(TMx(end-1)));
  title(txt,'Fontsize',11);
  
  set(gca,'tickdir','out',...
	  'xlim',[round(TMx(1))-0.01 TMx(end)+0.1],...
	  'xtick',[round(TMx(1)):round(TMx(end))],...
	  'ylim',[-1000 0]);
  
% Plot dlt S  
  axes('Position',[0.08 0.08 0.85 0.4]);
  pcolor(TMx,ZZ,dltS); shading flat
  caxis([cs1 cs2]);
  colormap(cmpS);
  CB = colorbar;
  set(CB,'Position',[0.935 0.08 0.01 0.4],...
	 'TickLength',0.04);
  
  txt = sprintf('%s, S',nm);
  title(txt,'Fontsize',11);
  
  set(gca,'tickdir','out',...
	  'xlim',[round(TMx(1))-0.01 TMx(end)+0.1],...
	  'xtick',[round(TMx(1)):round(TMx(end))],...
	  'ylim',[-1000 0]);
  
  bottom_text(txtb,'pwd',1);
  
  if s_fig>0
    fgnm=sprintf('%sdltTSvsTime_EN4_%s',pthfig,nm);
    fprintf('Saving %s\n',fgnm);
    print('-dpng','-r250',fgnm);
  end

end