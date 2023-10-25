% 0.04 HYCOM-CICE
%
% Calculate/extract T,S,U 
% Daily values
% across SE section on the Gr Shelf
% WOCE transect - see Sutherland and Pickart, PiO, 2008
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

dday=7; 
YR1=2019;
YR2=2019;
YRS = [YR1:YR2];
s_mat=2; 


ixx    = 9; % experiment name and dir - check with EXPT - expt 023
%ixx    = 6;  % expt 022 original 
EXPT   = sub_cice_experiments;
expt   = EXPT(ixx).Nmb;
texpt  = EXPT(ixx).cice_opt; % CICE options for sens. experiments
res    = EXPT(ixx).res;

Cp = 4200; % J/kg K
Tref1 = -1.8; % Ref T to calc. H flux
Tref2 = 0;    % Ref T to calc. H flux
Sref1 = 34.8;
Sref2 = 34.9;
rhow=1027; 
hgg=1e20;
f_zgrd=0;  % =1 - calculate fluxes from z-grid interpolated U,T,S - less accurate
           % mostly for comparison and validation

pthout = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.04/data_straits/';
%pthmat = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.04/data_theresa/';
pthmat = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/%3.3i/strait_fluxes/',expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';

regn = 'ARCc0.04';
fprintf('%s-%3.3i Daily Greenl Shelf SE Section U,T,S Zlevels, %i-%i\n\n',...
	regn,expt,YR1,YR2);

ftopo = sprintf('%s/depth_%s_17DD.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);

GC = sub_greenl_isobath_arc04(HH,LON,LAT);
Ig=GC.cntr_Iindx;
Jg=GC.cntr_Jindx;


% Specify sections:
% Greenlaaand    shelf OSNAP mooring line
SCT.Name='EastGrShelf_OSNAP';
SCT.IJ=[ 1226, 850; 1342, 810];
nsct=1;

fprintf('Section: %s\n',SCT.Name);
IJs=SCT.IJ;
[IIs,JJs]=sub_xsct_indx(IJs(1,1),IJs(1,2),IJs(2,1),IJs(2,2));
SCT.I=IIs;
SCT.J=JJs;

% Find normal at every small segments (+ northward, + eastward)
% Find intersection with GR contour
for ip=1:nsct
  IIs=SCT(ip).I;
  JJs=SCT(ip).J;
  ni=length(IIs);
  ist=IIs(1);
  jst=JJs(1);
  
  clear Nrm Hbtm dL Dst
  for isg=1:ni-1
    i0 = IIs(isg);
    i1 = IIs(isg+1);
    j0 = JJs(isg);
    j1 = JJs(isg+1);

    if i1==i0 % Y section, X norm ->
      Nrm(isg,1)=1;
      Nrm(isg,2)=0;  
    else    % Y norm ^
      Nrm(isg,1)=0;
      Nrm(isg,2)=1;
    end
    Hbtm(isg,1)=HH(j0,i0);
% Segments length:    
    dl=distance_spheric_coord(LAT(j0,i0),LON(j0,i0),LAT(j1,i1),LON(j1,i1));
    dL(isg,1)=dl;
    dst=distance_spheric_coord(LAT(j0,i0),LON(j0,i0),...
			       LAT(jst,ist),LON(jst,ist));
    Dst(isg,1)=dst;
  end
  Hbtm(ni)=HH(j1,i1);
  Dst(ni)=distance_spheric_coord(LAT(j1,i1),LON(j1,i1),...
			       LAT(jst,ist),LON(jst,ist));
  dL(ni)=dl;
  
  SCT(ip).Nrm=Nrm;
  SCT(ip).Hbottom=Hbtm;
% Intersection pnt:
  dd=sqrt((IIs-Ig).^2+(JJs-Jg).^2);
  [jm,im]=find(dd==min(min(dd)),1);
  imn=IIs(im);
  jmn=JJs(im);
  
  SCT(ip).GrCntr_intrcp=im;
  SCT(ip).GrCntr_I=imn;
  SCT(ip).GrCntr_J=jmn;
  SCT(ip).Segm_dL=dL;
  SCT(ip).Dist_origin=Dst;
  
end

btx = 'flux_SEGrShelf_xsct_gofs35_arc04.m';
f_map=0;
if f_map==1
  fn=1;
  sub_plot_Greenl_contour004(HH,LON,LAT,fn,GC);

  for ip=1:nsct
%    IJ=SCT(ip).IJ;
%    plot([IJ(1,1) IJ(2,1)],[IJ(1,2) IJ(2,2)],...
%	 'Linewidth',2.5,'Color',[1. 0.6 0]);
    IIs=SCT(ip).I;
    JJs=SCT(ip).J;
    plot(IIs,JJs,'-',...
	 'Linewidth',2.5,'Color',[1. 0.6 0]);
  end
  bottom_text(btx,'pwd',1);
end

% Interpolate onto fixed z-grid
% Specify fixed ZZ (interf) and ZM levels:
ZZf = [(0:-1:-20)';(-22:-2:-50)';(-55:-5:-500)';...
       (-510:-10:-2000)';(-2025:-25:-3000)';(-3050:-50:-5000)'];
kzz = length(ZZf);

dZf=diff(ZZf);
ZMf = [];
for ik=1:kzz-1
  ZMf(ik,1)=ZZf(ik)+0.5*(ZZf(ik+1)-ZZf(ik));
end

for YR=YR1:YR2
  yr=YR;
  fmatu=sprintf('%s%3.3i_SEGreenlSh_xsct_dayUTSZ_%i.mat',pthmat,expt,YR);
  
  UTSZ=struct;
  UTSZ = SCT;
  UTSZ.Title = 'U norm, m/s, SE Greenl shelf sections, Z levels, monthly mean';
  UTSZ.ZZlevels = ZZf;
  UTSZ.ZM = ZMf;

  cc=0; 
  for iday=1:dday:365
    dj1=datenum(yr,1,1);
    dnmb = dj1+iday-1; 
    DV   = datevec(dnmb);

    pthbin=sprintf('/nexsan/people/ddmitry/hycom/ARCc0.04_%3.3i/data/%4.4i_%s/',...
                    expt,yr,texpt);

    fina = sprintf('%s%3.3i_archv.%4.4i_%3.3i_00.a',pthbin,expt,yr,iday);
    finb = sprintf('%s%3.3i_archv.%4.4i_%3.3i_00.b',pthbin,expt,yr,iday);

    if yr==2017 & DV(2)<4
      pthbin='/nexsan/people/ddmitry/hycom/ARCc0.04_022/data/2017_BL99Tfrz/';
      fina = sprintf('%s022_archv.%4.4i_%3.3i_00.a',pthbin,yr,iday);
      finb = sprintf('%s022_archv.%4.4i_%3.3i_00.b',pthbin,yr,iday);
    end

    if ~exist(fina,'file') | ~exist(finb,'file')
      fprintf('Not found %s or %s\n\n',fina,finb);
      continue;
    end

    cc=cc+1;

    fprintf('Reading %4.4i/%2.2i/%2.2i: %s\n',DV(1:3),fina);

    tic;
    [F,n,m,nlr] = read_hycom(fina,finb,'temp');
    F(F>hgg)=nan;
    T=F;

    [F,n,m,nlr] = read_hycom(fina,finb,'salin');
    F(F>hgg)=nan;
    S=F;
%
    [F,n,m,nlr] = read_hycom(fina,finb,'u-vel.');
    F(F>hgg)=0;
    U=F;

    [F,n,m,llr] = read_hycom(fina,finb,'u_btrop');
    F(F>hgg)=0;
    Ub=squeeze(F);

    for ilr=1:nlr
      U(ilr,:,:)=squeeze(U(ilr,:,:))+Ub;
    end

% These are NOT total velocities    - archv
% barotropic needed
    [F,n,m,nlr] = read_hycom(fina,finb,'v-vel.');
    F(F>hgg)=0;
    V=F;

    [F,n,m,llr] = read_hycom(fina,finb,'v_btrop');
    F(F>hgg)=0;
    Vb=squeeze(F);

    for ilr=1:nlr
      V(ilr,:,:)=squeeze(V(ilr,:,:))+Vb;
    end

    [ZMh,ZZh] = sub_zz_zm(fina,finb,HH);
    dH=abs(diff(ZZh,1));

    AA = sub_sct_interp2z(SCT,U,V,T,S,ZZf,ZZh,HH,dH,DX,DY);    
%keyboard
    UTSZ.Unrm(cc,:,:) = AA.NormalU;
    UTSZ.Temp(cc,:,:) = AA.Temp;
    UTSZ.Saln(cc,:,:) = AA.Salin;
    UTSZ.Time(cc)     = dnmb;	

    fprintf('++++++>  Processed 1 record %8.5f min\n\n',toc/60);

    if mod(cc,5)==0 & s_mat>0
      fprintf('Saving %s\n',fmatu);
      save(fmatu,'UTSZ');
    end	

  end  % dday
  
  if s_mat>0
    fprintf('Saving %s\n',fmatu);
    save(fmatu,'UTSZ');
  end	
  
end;   % years    
    
%if s_mat>0
%  fprintf('Saving %s\n',fmatu);
%  save(fmatu,'UTSZ');
%end	
  
  
  
  
