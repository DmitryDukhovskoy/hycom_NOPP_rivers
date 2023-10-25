% To test hypothesis about convergence/divergence
% of particles in the SPNA at different depths
%
% Estimate divergence - volume averaged
% using Gauss (diverg) theorem: integrate fluxes over the 
% boundaries
%
%
addpath /Net/Movies0/ddmitry/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /home/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers;
%addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
%addpath /Net/Movies0/ddmitry/MyMatlab/seawater
startup

close all
clear

s_mat = 2012; % = 1 - save from beginning
              % = 2003 - upload saved mat file 2003, find last saved record and start from next step
              % = 0 - do not save

expt = 112;
TV   = 11;
YR1  = 2012;
YR2  = 2012;

hgg = 1e25;

pthtopo= '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/110/data_GG_prt/';

btx='calc_flux_SPNA.m';

fprintf('arc08-%3.3i UV transport SPNA %i-%i, save=%i\n',...
 expt,YR1,YR2,s_mat);


ftopo = sprintf('%sdepth_ARCc0.08_%2.2i.nc',pthtopo,TV); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);

load('SPG_noNorth_indx.mat');
IGR(end+1,:)=IGR(1,:);

%
% Get indices of segments
clear SCT
nsct=length(IGR)-1;
SCT.IGR=IGR;
for ii=1:nsct
  SCT(ii).IJ = [IGR(ii,:); IGR(ii+1,:)];
  IJs=SCT(ii).IJ;
  [IIs,JJs]=sub_xsct_indx(IJs(1,1),IJs(1,2),IJs(2,1),IJs(2,2));
  SCT(ii).I=IIs;
  SCT(ii).J=JJs;
  nsg=length(IIs);
  clear XX YY Hb
  for ij=1:nsg
    i0=IIs(ij);
    j0=JJs(ij);
    XX(ij)=LON(j0,i0);
    YY(ij)=LAT(j0,i0);
    Hb(ij)=HH(j0,i0);
  end
  SCT(ii).long=XX;
  SCT(ii).latd=YY;
  SCT(ii).Hb=Hb;
  iocn = length(find(Hb<0));
  SCT(ii).Nocn = iocn; % # of ocean points in the segment
end


f_map=1;
if f_map==1
  figure(10); clf;
  hold on;

  lcmp=[0 0 0; 1 1 1];
  Lmsk=HH*0;
  Lmsk(HH<0)=1;
  pcolor(Lmsk); shading flat;
  colormap(lcmp);

  contour(HH,[-5000:500:0],'Color',[0.8 0.8 0.8]);
  caxis([0 1]);

  for isg=1:nsct
    II = SCT(isg).I;
    JJ = SCT(isg).J;

    iocn = SCT(isg).Nocn;
    if iocn<=0
      clr = [0.7 0.7 0.7];
    else
      clr = [1 0.4 0];
    end

    II = SCT(isg).I;
    JJ = SCT(isg).J;

    plot(II,JJ,'.-','Color',clr);
  end

  axis('equal');
  set(gca,'xlim',[300 1200],...
            'ylim',[80 1100]);

  title('SPNA region');
  bottom_text(btx,'pwd',1);

  fsct = sprintf('%sSCT_SPNA.mat',pthmat);
  fprintf('Saving Sections SPNA %s\n',fsct);
  save(fsct,'SCT');

keyboard
end

% Find contour through UV pnts for SPNA contour
% Not collocated U&V
Ig=[];
Jg=[];
for isc=1:nsct
  II = SCT(isc).I;
  JJ = SCT(isc).J;
  II = II(:);
  JJ = JJ(:);
  if isc>1
    dd=sqrt((Ig(end)-II(1))^2+(Jg(end)-JJ(1)));
    if dd==0, II=II(2:end); JJ=JJ(2:end); end;
  end
  Ig = [Ig;II];
  Jg = [Jg;JJ];
end
nin  = 1; % positive inside the contour
UVGR = sub_UVpnts_contour(Ig,Jg,nin,HH);
dL   = sub_segm_dL(UVGR.sgmX,UVGR.sgmY,LON,LAT);

GC.Norm=UVGR.Norm;
GC.sgmX=UVGR.sgmX;
GC.sgmY=UVGR.sgmY;
GC.Ictr=UVGR.Ictr;
GC.IJ_indx=UVGR.gridIndx_IJ; % grid indices corresponding to UV segments
GC.IJ_adj_indx=UVGR.adjIndx_I1J1; % adjacent grid point index - for collocating dH to UV pnts
GC.segm_dL=dL;

f_sct=0;
if f_sct==1
  nf=15;
  sub_check_sct(nf,HH,GC);
end


for YR=YR1:YR2
  yr=YR;
  dE=datenum(yr,12,31);
  dJ1=datenum(yr,1,1);
  ndays=dE-dJ1+1;

  cc=0;


% Daily mean vol fluxes
  GVHFLX = struct;
  dday  = 7;
  tm0   = 0;
  TM    = [];
  

  if YR<s_mat
    fprintf('Starting from saved year %i, skipping %i\n',s_mat,YR);
    continue
  elseif YR==s_mat
    fmatout=sprintf('%shycom008_%3.3i_SPNAcntr_VolFlux_%4.4i.mat',...
                    pthmat,expt,YR);
    fprintf('Loading %s\n',fmatout);
    load(fmatout);

    cc  = length(GVHFLX.Time);
    tm0 = GVHFLX.Time(cc);

    fprintf('Last saved record: cc=%i, %s\n',cc,datestr(tm0));
    fprintf('Continue from %s\n\n',datestr(tm0+dday));

  end


  for iday=1:dday:ndays;
    dnmb=datenum(dJ1)+iday-1;
    DV   = datevec(dnmb);
    iday = dnmb-datenum(yr,1,1)+1;

    if dnmb<=tm0, 
      fprintf('RESTART: Skipping %s\n',datestr(dnmb));
      continue; 
    end;

%    keyboard

    pthbin = sprintf('/nexsan/archive/ARCc0.08_%3.3i/data/%4.4i/',expt,yr);
    if expt==112,
      pthbin=sprintf('/nexsan/hycom/ARCc0.08_%3.3i/data/%4.4i/',expt,yr);
    end

    fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
    finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);

    if ~exist(fina,'file') | ~exist(finb,'file')
      fprintf('Not found %s or %s\n\n',fina,finb);
      continue;
    end

    cc=cc+1;

    fprintf('Reading %4.4i/%2.2i/%2.2i: %s\n',DV(1:3),fina);

    tic;
%    [F,n,m,nlr] = read_hycom(fina,finb,'temp');
%    F(F>hgg)=nan;
%    T=F;

%    [F,n,m,nlr] = read_hycom(fina,finb,'salin');
%    F(F>hgg)=nan;
%    S=F;
%
    [F,n,m,nlr] = read_hycom(fina,finb,'u-vel.');
    F(F>hgg)=0;
    U=F;

% These are total velocities    - mean output fields
% barotropic not needed
%
%    [F,n,m,l] = read_hycom(fina,finb,'u_btrop');
%    F(F>hgg)=nan;
%    u_btrop=squeeze(F);

    [F,n,m,nlr] = read_hycom(fina,finb,'v-vel.');
    F(F>hgg)=0;
    V=F;

    [ZMh,ZZh] = sub_zz_zm(fina,finb,HH);
    dH=abs(diff(ZZh,1));

%  =========================================
% Calculate fluxes through SPNA Contour
%  =========================================
    iGRs=GC.sgmX;
    jGRs=GC.sgmY;
    nns=length(iGRs(:,1));
    fprintf('Fluxes through SPNA Contour ...\n');
    for j=1:nns
      if mod(j,400)==0, fprintf('  done %6.2f%%...\n',j/nns*100); end
      x1=iGRs(j,1);
      x2=iGRs(j,2);
      y1=jGRs(j,1);
      y2=jGRs(j,2);
      xnrm=GC.Norm(j,1);
      ynrm=GC.Norm(j,2);
% Current/adj cells - collocation
      i0=GC.IJ_indx(j,1);
      j0=GC.IJ_indx(j,2);
      i1=GC.IJ_adj_indx(j,1);
      j1=GC.IJ_adj_indx(j,2);

      if HH(j0,i0)>=0 & HH(j1,i1)>=0
        GUTS.Unrm(:,j) = zeros(nlr,1);
        GUTS.dH(:,j)   = zeros(nlr,1);
        GUTS.Hb(j)   = 0;
        continue;
      end

% Collocate H to U-V points
      T=[];
      S=[];
      if y1==y2  % horiz segment, V flux
        CLC=sub_collocate_H2uv(V,dH,HH,i0,j0,i1,j1);
      else
        CLC=sub_collocate_H2uv(U,dH,HH,i0,j0,i1,j1);
      end
      In=find(CLC.dHn<1e-3);
%      CLC.Tn(In)=nan;
%      CLC.Sn(In)=nan;
      CLC.Un(In)=nan;

      snrm=sign(xnrm+ynrm); % one is 0,  need only direction
      GUTS.Hb(j)=-CLC.Hn;
      GUTS.Unrm(:,j)=snrm*CLC.Un;
%      GUTS.Tnrm(:,j)=CLC.Tn;
%      GUTS.Snrm(:,j)=CLC.Sn;
      GUTS.dH(:,j)=CLC.dHn;

      chgr=0;
      if chgr==1
        plot([x1 x2],[y1 y2],'b-');
        plot(i0,j0,'r*');
        plot(i1,j1,'gd');
        xp=0.5*(x1+x2);
        yp=0.5*(y1+y2);
        plot([xp xp+xnrm],[yp yp+ynrm],'g-');
      end
    end
    un=GUTS.Unrm;
%    tn=GUTS.Tnrm;
%    sn=GUTS.Snrm;
    DZ=abs(GUTS.dH);
    dx=GC.segm_dL;
    [DX,dmm]=meshgrid(dx,[1:nlr]);
%    rhow=sw_dens0(sn,tn);
%keyboard
% Fluxes over all sections    
    vfG0=nansum(un.*DX.*DZ); % total trunsport
    GVHFLX.Time(cc)=dnmb;
    GVHFLX.Vflx(cc)=nansum(vfG0); % overall vol flux
    
%
%   Fluxes pointwise across the SPNA contour:
    dmm = un.*DX.*DZ;  % transport in grid cells m3/s
    GVHFLX.VflxPnts(cc,:,:)=dmm;
    GVHFLX.DZ(cc,:,:)=DZ;
%
%  Fluxes averaged over the layers:
% 10, 15, 23, 31
    vf1 = nansum(dmm(1:10,:)); 
    vf2 = nansum(dmm(11:15,:));
    vf3 = nansum(dmm(16:23,:));
    vf4 = nansum(dmm(24:31,:));

% Mean:
    if cc==1
      GVHFLX.Vfm = 0;
      GVHFLX.Vf1 = 0; 
      GVHFLX.Vf2 = 0; 
      GVHFLX.Vf3 = 0; 
      GVHFLX.Vf4 = 0; 
    end
    GVHFLX.Vfm = GVHFLX.Vfm+nansum(vfG0);
    GVHFLX.Vf1 = GVHFLX.Vf1+nansum(vf1); % 0-50
    GVHFLX.Vf2 = GVHFLX.Vf2+nansum(vf2); % lr 50-150
    GVHFLX.Vf3 = GVHFLX.Vf3+nansum(vf3); % lr 50-150
    GVHFLX.Vf4 = GVHFLX.Vf4+nansum(vf4); % lr 50-150

% Mean fluxes:
    Vfm0  = GVHFLX.Vfm/cc;
    Vfm1  = GVHFLX.Vf1/cc;
    Vfm2  = GVHFLX.Vf2/cc;
    Vfm3  = GVHFLX.Vf3/cc;
    Vfm4  = GVHFLX.Vf4/cc;


    fprintf('Mean Vol Fluxes (+ into box)\n');
    fprintf(' TOTAL     VOL = %5.1f Sv, \n',Vfm0*1e-6);
    fprintf(' Lr 01:10  Vol = %5.1f Sv, \n',Vfm1*1e-6);
    fprintf(' Lr 11:15  Vol = %5.1f Sv, \n',Vfm2*1e-6);
    fprintf(' Lr 16:23  Vol = %5.1f Sv, \n',Vfm3*1e-6);
    fprintf(' Lr 24:31  Vol = %5.1f Sv, \n',Vfm4*1e-6);

    fprintf('++++++>  Processed 1 record %8.5f min\n\n',toc/60);


    fmatout=sprintf('%shycom008_%3.3i_SPNAcntr_VolFlux_%4.4i.mat',...
                     pthmat,expt,YR);

    if s_mat>0 & mod(cc,10)==0
      fprintf('Saving %s\n',fmatout);
      save(fmatout,'GVHFLX');
    end


  end
  if s_mat>0
		  fprintf('Saving %s\n',fmatout);
		  save(fmatout,'GVHFLX');
  end

end

%if s_mat>0
%		fprintf('Saving %s\n',fmatout);
%		save(fmatout,'GVHFLX');
%end



