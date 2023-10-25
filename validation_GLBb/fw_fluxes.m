% Calculate FW fluxes - annual mean
% across straits 
% /nexsan/GLBa0.08/expt_90.9/data
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

s_mat  = 1; % 
s_fig  = 1;


rg=9806;  % convert pressure to depth, m
%Sref=34.8;
Sref=35;

pth190='/nexsan/GLBb0.08/GLBb0.08_190/data/meanstd/'; % 1993-1994
pth191='/nexsan/GLBb0.08/GLBb0.08_191/data/meanstd/'; % 1995-2012
pthmat = '/Net/mars/ddmitry/hycom/GLBb0.08/data_mat/';
ftopo = sprintf('%sGLBb_T07_subset_Natl.mat',pthmat);
pthglb = '/nexsan/GLBb0.08/GLBb0.08_191/topo/';
fmat = sprintf('%sFWflx_annual_Sref%3i_1993-2015.mat',pthmat,round(Sref*10));

YRPLT=[];
cnc=0;
for yr=1993:2015
  cnc=cnc+1;
  YRPLT(cnc,1)=yr;
end

fprintf(' Save mat =%i, %s\n',s_mat,fmat);
fprintf('  Time: %i - %i\n',YRPLT(1),YRPLT(end));

f_getH=0; % prepare sub-region
if f_getH>0
  expt=190;
  sub_get_GLBgrid(pthglb,expt,ftopo);
end
fprintf('Loading topo %s\n',ftopo);
load(ftopo);

% Sections (2 pnts / section):
IJs=[   145,   789;
   298,   789;
   151,   525;
   247,   525;
    49,   302;
   238, 302;
   604,  807;
   749,   807;
   383,   397;
   532,   397;
   621, 357;
   889,   357;
   121, 11;
   819, 11];
nvv=size(IJs,1);
NMS{1}='BaffinBay';
NMS{2}='DavisStr';
NMS{3}='LabrSea';
NMS{4}='FramStr';
NMS{5}='DenmarkStr';
NMS{6}='Fareo';
NMS{7}='NAtl';

cc=0;
for ii=1:2:nvv-1
  cc=cc+1;
  SEGM(cc).IJ(1,1:2)=IJs(ii,1:2);
  SEGM(cc).IJ(2,1:2)=IJs(ii+1,1:2);
end;

fprintf('# of segments: %i\n',cc);

p_segm=0;
if p_segm==1
  figure(10);
  clf;
  contour(HH,[0 0],'k');
  hold on
  for ii=1:2:nvv-1
    plot([IJs(ii,1) IJs(ii+1,1)],[IJs(ii,2) IJs(ii+1,2)],'r');
  end
end


% Mask of the region of interest
% exclude Pacific Ocean and North.Atl.
LMSK=HH;
LMSK(LMSK>=0)=0;
LMSK(LMSK<0)=1;
IDEEP=find(LMSK>0);
INAN =find(LMSK==0);

% Grid cell spacing
[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2

% Create long-term means from yearly means:
%for ik=1:2
np=length(YRPLT);
iy1=1;
cc=iy1-1;
for iyr=iy1:np
  tic;
  yr=YRPLT(iyr);
  fprintf('Calculating FWC for yr=%i\n',yr);
  
  if yr<1995
    pthbin=pth190;
    expt=190;
  elseif yr>=1995 & yr <2013
    pthbin=pth191;
    expt=191;
  elseif yr==2013
    pthbin='/nexsan/GLBa0.08/expt_91.0/data/meanstd/';
    expt=910; % 
  elseif yr==2014 | yr==2015
    pthbin='/nexsan/GLBa0.08/expt_91.1/data/meanstd/';
    expt=911; % 
  end
  
  fnm=sprintf('%s%3.3i_archMN.%4.4i_01_%4.4i_12',pthbin,expt,yr,yr);
  fina=[fnm,'.a'];
  finb=[fnm,'.b'];
  
  ind1=IND.i1;
  ind2=IND.i2;
  jnd1=IND.j1;
  jnd2=IND.j2;
  
  fld='thknss';
  [F,n,m,l] = read_hycom(fina,finb,fld);
  F(F>1e20)=nan;
  DP=F(:,jnd1:jnd2,ind1:ind2)./rg; % subset
  DP(DP<0.1)=nan; % 0-m layers, vanished

% Interf. depths and depths of the middle of the layers (m)
% NOTE: sign convention: depths are negative
% dP - Pa or m
  [ZZ,ZM] = sub_thck2dpth(DP); % 

  fld='salin';
  [F,n,m,l] = read_hycom(fina,finb,fld);
  F(F>1e20)=nan;
  S=F(:,jnd1:jnd2,ind1:ind2);

%  fld='u-vel.';
%  [F,n,m,l] = read_hycom(fina,finb,fld);
%  F(F>1e20)=nan;
%  U=F(:,jnd1:jnd2,ind1:ind2);

  fld='v-vel.';
  [F,n,m,l] = read_hycom(fina,finb,fld);
  F(F>1e20)=nan;
  V=F(:,jnd1:jnd2,ind1:ind2);

  f_btrop=0;
  if f_btrop>0
    fld='u_btrop';
    [F,n,m,l] = read_hycom(fina,finb,fld);
    F(F>1e20)=nan;
    F=squeeze(F);
    Ub=F(:,jnd1:jnd2,ind1:ind2);

    fld='v_btrop';
    [F,n,m,l] = read_hycom(fina,finb,fld);
    F(F>1e20)=nan;
    F=squeeze(F);
    Vb=F(jnd1:jnd2,ind1:ind2);
  end

%  keyboard
  
  f_chck=0;
  if f_chck>0, % 
    i0=500;
    j0=100;
    u=squeeze(U(:,j0,i0));
    v=squeeze(V(:,j0,i0));
    dp=squeeze(DP(:,j0,i0));
    h=nansum(dp);
    uav=nansum(u.*dp)./h;
    vav=nansum(v.*dp)./h;
    ub=Ub(j0,i0);
    vb=Vb(j0,i0);
    if abs(1-uav/ub)<0.1 & abs(1-vav/vb)<0.1
      fprintf('u-vel, v-vel is total\n');
    else 
      fprintf('u-vel is baroclinnnic\n');
    end
    
  end;
 
  cc=cc+1;
  S(:,INAN)=nan;
  NS=length(SEGM);
  for isgm=1:NS
    i1=SEGM(isgm).IJ(1,1);
    i2=SEGM(isgm).IJ(2,1);
    j1=SEGM(isgm).IJ(1,2);
    j2=SEGM(isgm).IJ(2,2);
    
    dp=squeeze(DP(:,j1:j2,i1:i2));
    ss=squeeze(S(:,j1:j2,i1:i2));
    if j1==j2
      uu=squeeze(V(:,j1:j2,i1:i2));
      dx=DX(j1:j2,i1:i2);
    else
      uu=squeeze(U(:,j1:j2,i1:i2));
      dx=DY(j1:j2,i1:i2)';
    end
    zm=squeeze(ZM(:,j1:j2,i1:i2));
% Follow ~Haine et al., 2015 definition of FWC
    [ll,pp]=size(ss);
    FWflx=zeros(ll,pp);
%keyboard
%  tic;
    for k=1:ll-1
      dz=abs(dp(k,:));
      dz(isnan(dz))=0;
      ss1=ss(k,:);
      ss2=ss(k+1,:);
      uu1=uu(k,:);
      uu2=uu(k+1,:);
      Ib1=find(ss1>Sref);
      ss1(Ib1)=nan;
      ss2(Ib1)=nan;
      Ib2=find(ss2>Sref); % CHeck next layer see if S > Sref there

      fwc=dz.*(Sref-ss1)/Sref; % m
      fwf=fwc.*dx.*uu1; %m3/s
% Take care about layers where S becomes S>Sref
% add FWC below layer interface from depth of Sref
% up to the bottom lyaer k, S=0.5(S(k)+Sref)
%    if ~isempty(Ib2)
%      zm1 = squeeze(ZM(k,j1:j2,i1:i2));
%      zm2 = squeeze(ZM(k+1,j1:j2,i1:i2));
%      dZm = zm2-zm1;
%      dS = ss2-ss1;
%      zmSref = zm1+dZm./dS.*(Sref-ss1);
%      zz = squeeze(ZZ(k+1,j1:j2,i1:i2)); % lower interface
%      dltZs  = abs(zz-zmSref);
%      Smn=0.5*(ss1+Sref); % for integration from 
%      dltFwc = dltZs.*(Sref-Smn)/Sref;
%      fwc(Ib2) = dltFwc(Ib2);
%    end
      FWflx(k,:)=fwf;
    end % layers
   [nl, nx]=size(uu);
   [X,Y]=meshgrid([1:nx],[1:nl]); 
%    keyboard
    
    FWFLX(cc).Title = 'Annual FW flux thr segments  from GLBb0.08';
    FWFLX(cc).Source=pthbin;
    FWFLX(cc).Year  = yr;
    FWFLX(cc).SEGM(isgm).Name=NMS{isgm};
    FWFLX(cc).SEGM(isgm).IJ=SEGM(isgm).IJ;
    FWFLX(cc).SEGM(isgm).FWFlux_m3 = FWflx;
    FWFLX(cc).SEGM(isgm).S_section=ss;
    FWFLX(cc).SEGM(isgm).X_section=X;
    FWFLX(cc).SEGM(isgm).ZM_section=zm;
    FWFLX(cc).Sref = Sref;
  end  % segments
%  Fwc(Ilnd)=nan; % exclude not needed regions
  fprintf('FWC calculation %6.2f min\n',toc/60);
  
%keyboard
% saving data
  if s_mat>0
    fprintf('Saving %s\n\n',fmat);
    save(fmat,'FWFLX');
  end;
end; % year



