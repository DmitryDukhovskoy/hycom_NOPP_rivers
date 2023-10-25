% Calculate heat & salt fluxes
% across sections on the Gr Shelf
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

YR1=2012;
YR2=2012;
YRS = [YR1:YR2];
s_mat=1; 


rg = 9806;
rhow=1027; 
hgg=1e20;

regn = 'ARCc0.08';
%expt = 110;
expt = 112; % Greenland runoff
pthfig  = '/nexsan/people/ddmitry/Net_ocean/hycom/ARCc0.08/110/fig_meanUV/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat =sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
pthmat2 =sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_GrSect/',expt);
btx='flux_Grshelf_sections.m';

fprintf('%s-%3.3i Monthly Greenl Shelf Sections U,T,S Zlevels, %i-%i\n\n',...
	regn,expt,YR1,YR2);

ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);

GC = sub_greenl_isobath(HH,LON,LAT);
Ig=GC.cntr_Iindx;
Jg=GC.cntr_Jindx;

% Specify sections:
% Greenlaaand    shelf OSNAP mooring line
SCT(1).Name='DenmarkStr';
SCT(1).IJ=[788 610; 847 536];
SCT(2).Name='EastGrShelf_OSNAP';
SCT(2).IJ=[   613   425; 671   405];
SCT(3).Name='DavisStr';
SCT(3).IJ=[471 673; 556 652];
SCT(4).Name='NaresStr';
SCT(4).IJ=[570 1037; 605 999];
nsct=length(SCT);

for ip=1:nsct
  fprintf('Section: %s\n',SCT(ip).Name);
  IJs=SCT(ip).IJ;
  [IIs,JJs]=sub_xsct_indx(IJs(1,1),IJs(1,2),IJs(2,1),IJs(2,2));
  SCT(ip).I=IIs;
  SCT(ip).J=JJs;
end

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


f_map=0;
if f_map==1
  fn=1;
  sub_plot_Greenl_contour(HH,LON,LAT,fn,GC);

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
ZZf = [(0:-1:-10)';(-12:-2:-30)';(-35:-5:-100)';...
       (-110:-10:-1000)';(-1025:-25:-2500)';(-2550:-50:-5000)'];
kzz = length(ZZf);

dZf=diff(ZZf);
ZMf = [];
for ik=1:kzz-1
  ZMf(ik,1)=ZZf(ik)+0.5*(ZZf(ik+1)-ZZf(ik));
end

UZGR = SCT;
UZGR(1).Title = 'U norm, m/s, Greenl shelf sections, Z levels, monthly mean';
UZGR(1).ZZlevels = ZZf;
UZGR(1).ZM = ZMf;

TZGR = SCT;
TZGR(1).Title = 'T, Greenl shelf sections, Z levels, monthly mean';
TZGR(1).ZZlevels = ZZf;
TZGR(1).ZM = ZMf;

SZGR = SCT;
SZGR(1).Title = 'S, Greenl shelf sections, Z levels, monthly mean';
SZGR(1).ZZlevels = ZZf;
SZGR(1).ZM = ZMf;

UVZGR = SCT;
UVZGR(1).Title = 'U,V components, Greenl shelf sections, Z levels, monthly mean';
UVZGR(1).ZZlevels = ZZf;
UVZGR(1).ZM = ZMf;


dday=7;
for YR=YR1:YR2
  yr=YR;
  fmatu=sprintf('%sarc08_expt%3.3i_greenl_shsect_mnthUzlv_%i.mat',pthmat,expt,YR);
  fmatt=sprintf('%sarc08_expt%3.3i_greenl_shsect_mnthTzlv_%i.mat',pthmat,expt,YR);
  fmats=sprintf('%sarc08_expt%3.3i_greenl_shsect_mnthSzlv_%i.mat',pthmat,expt,YR);
  fmatuv=sprintf('%sarc08_expt%3.3i_greenl_shsect_mnthUVz_%i.mat',pthmat,expt,YR);

  for imo=1:12
    d1=datenum(yr,imo,1);
    d2=d1+32;
    dv2=datevec(d2);
    dm=datenum(dv2(1),dv2(2),1)-d1;

    for ip=1:nsct
      I=SCT(ip).I;
      np=length(I);
      TMP(ip).Usm = zeros(kzz-1,np);
      TMP(ip).Tsm = zeros(kzz-1,np);
      TMP(ip).Ssm = zeros(kzz-1,np);
      TMP(ip).UUsm= zeros(kzz-1,np);
      TMP(ip).VVsm= zeros(kzz-1,np);
    end

    cc=0;
    for mdd=1:dday:dm
      dnmb = datenum(yr,imo,mdd);
      DV   = datevec(dnmb);
      iday = dnmb-datenum(yr,1,1)+1;
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

      [F,n,m,nlr] = read_hycom(fina,finb,'v-vel.');
      F(F>hgg)=0;
      V=F;


      [ZMh,ZZh] = sub_zz_zm(fina,finb,HH);
      dH=abs(diff(ZZh,1));
      UTSZ = sub_sct_interp2z(SCT,U,V,T,S,ZZf,ZZh,HH,dH,DX,DY);    

      for ip=1:nsct
	Usm=TMP(ip).Usm;
	Uz=UTSZ(ip).NormalU;
	TMP(ip).Usm = Usm+Uz;
	
	Tsm=TMP(ip).Tsm;
	Tz=UTSZ(ip).Temp;
	TMP(ip).Tsm=Tsm+Tz;
	
        Ssm=TMP(ip).Ssm;
	Sz=UTSZ(ip).Salin;
	TMP(ip).Ssm = Ssm+Sz;
	
	usm=TMP(ip).UUsm;
	uz=UTSZ(ip).Ucomp;
	TMP(ip).UUsm=usm+uz;
	
	vsm=TMP(ip).VVsm;
	vz=UTSZ(ip).Vcomp;
	TMP(ip).VVsm=vsm+vz;
      end

      fprintf('++++++>  Processed 1 record %8.5f min\n\n',toc/60);

    end  % month day

    for ip=1:nsct
      UZGR(ip).Nrec(imo)=cc;
      Usm=TMP(ip).Usm;
      Usm = Usm/cc;
      UZGR(ip).U(imo,:,:)=Usm;
      
      TZGR(ip).Nrec(imo)=cc;
      Tsm=TMP(ip).Tsm;
      Tsm=Tsm/cc;
      TZGR(ip).T(imo,:,:)=Tsm;
      
      SZGR(ip).Nrec(imo)=cc;
      Ssm=TMP(ip).Ssm;
      Ssm=Ssm/cc;
      SZGR(ip).S(imo,:,:)=Ssm;
      
      UVZGR(ip).Nrec(imo)=cc;
      usm=TMP(ip).UUsm;
      usm=usm/cc;
      UVZGR(ip).Ucomp(imo,:,:)=usm;
      
      vsm=TMP(ip).VVsm;
      vsm=vsm/cc;
      UVZGR(ip).Vcomp(imo,:,:)=vsm;
      
      umx = max(max(Usm));
      umn = min(min(Usm));
      tmx = max(max(Tsm));
      tmn = min(min(Tsm));
      smx = max(max(Ssm));
      smn = min(min(Ssm));
      
      fprintf('Section: %s\n',SCT(ip).Name);
      fprintf('End of Month: umx = %6.3f umn=%6.3f\n',umx,umn);
      fprintf('End of Month: tmx = %6.3f tmn=%6.3f\n',tmx,tmn);
      fprintf('End of Month: smx = %6.3f smn=%6.3f\n',smx,smn);

    end  % dday
    
    if s_mat==1
      fprintf('Saving %s\n',fmatu);
      save(fmatu,'UZGR');

      fprintf('Saving %s\n',fmatt);
      save(fmatt,'TZGR');

      fprintf('Saving %s\n',fmats);
      save(fmats,'SZGR');
      
      fprintf('Saving %s\n',fmatuv);
      save(fmatuv,'UVZGR');
    end  

  end    % months
    
end;   % years    
    
  
  
  
  