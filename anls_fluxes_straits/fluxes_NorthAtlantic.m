% Calculate vol, heat, salt fluxes
% through N. Atlantic openings
% for budgets
%
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/Wavelet;
startup;

close all
clear


expt=112;
TV=11;
YR1=2007;
YR2=2007;

s_mat=1;

Cp = 4200; % J/kg K
Tref1= -1.8; % Ref T to calc. H flux
Sr1=34.8;
Sr2=34.95;
hgg=1e20;
f_zgrd=0;  % =1 - calculate fluxes from z-grid interpolated U,T,S - less accurate
           % mostly for comparison and validation

pthmat = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_natl/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';

btx='fluxes_NorthAtlantic.m';

fprintf('arc08-%3.3i Heat and Vol fluxes North Atlantic %i-%i, save=%i\n',...
	expt,YR1,YR2,s_mat);


ftopo = sprintf('%sdepth_ARCc0.08_%2.2i.nc',pthtopo,TV); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);


% Get gates coordinates from POP
SCT(1).Name='DavisStr';
SCT(1).IJ=[471 673; 556 652];
SCT(2).Name='Hudson';
SCT(2).IJ=[375, 523; 375, 577];
SCT(3).Name='NAtl';
SCT(3).IJ=[435, 251; 1068, 150];
SCT(4).Name='BarentsSea';
SCT(4).IJ=[1141,887; 1265 751];
SCT(5).Name='FramStr';
SCT(5).IJ=[917, 930; 1087, 930];
SCT(6).Name='NaresStr';
SCT(6).IJ=[570 1037; 605 999];
SCT(7).Name='LancasterSound';
SCT(7).IJ=[475 1017; 507 1017];
nsct=length(SCT);

%
% Orientation points to determine 
% +/- fluxes across the sections
% Poisitive is northward or towards the AO
IJPR=[525, 787;...
     525, 787;...
     784,392;...
     1294,987;...
     995, 1159;...
     995,1159;...
     492, 1296];
ni=length(IJPR);
for ii=1:ni
  i1=IJPR(ii,1);
  j1=IJPR(ii,2);
  IPR(ii)=sub2ind(size(HH),j1,i1);
end


for ip=1:nsct
  fprintf('Section: %s\n',SCT(ip).Name);
  IJs=SCT(ip).IJ;
  [IIs,JJs]=sub_xsct_indx(IJs(1,1),IJs(1,2),IJs(2,1),IJs(2,2));
  SCT(ip).I=IIs;
  SCT(ip).J=JJs;
  nsg=length(IIs);
  clear XX YY
  for ii=1:nsg
    i0=IIs(ii);
    j0=JJs(ii);
    XX(ii)=LON(j0,i0);
    YY(ii)=LAT(j0,i0);
  end
  SCT(ip).long=XX;
  SCT(ip).latd=YY;
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
  SCT(ip).Segm_dL=dL;
  SCT(ip).Dist_origin=Dst;
  
end


f_map=0;
if f_map==1
  fn=1;
%  sub_plot_Greenl_contour(HH,LON,LAT,fn,GC);
  figure(1); clf;
  contour(HH,[0 0],'k');
  hold on;
  contour(HH,[-5000:500:-100],'Color',[0.6 0.6 0.6]);
 
  for ip=1:nsct
%    IJ=SCT(ip).IJ;
%    plot([IJ(1,1) IJ(2,1)],[IJ(1,2) IJ(2,2)],...
%	 'Linewidth',2.5,'Color',[1. 0.6 0]);
    IIs=SCT(ip).I;
    JJs=SCT(ip).J;
    plot(IIs,JJs,'-',...
	 'Linewidth',2.5,'Color',[1. 0.3 0]);
  end
  
  axis('equal');
  set(gca,'xlim',[450 1100],...
	  'ylim',[350 1200]);
  
  bottom_text(btx,'pwd',1);
end


% Daily Heat, Salt, and vol fluxes
VHFLX=SCT;
for isc=1:nsct
  VHFLX(isc).Tref1=Tref1;
  VHFLX(isc).Sref1=Sr1;
  VHFLX(isc).Sref2=Sr2;
end

for YR=YR1:YR2
  yr=YR;
  dE=datenum(yr,12,31);
  dJ1=datenum(yr,1,1);
  ndays=dE-dJ1+1;
  
  cc=0;
  for isc=1:nsct
    VHFLX(isc).Vfm=0;  % track overal means, vol flux
    VHFLX(isc).Hfm1=0; % heat flx
    VHFLX(isc).Sfm1=0;
    VHFLX(isc).Sfm2=0;
  end
  
  for iday=1:ndays
    dnmb=datenum(dJ1)+iday-1;
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

% These are total velocities    
% barotropic not needed
%
%    [F,n,m,l] = read_hycom(fina,finb,'u_btrop');
%    F(F>hgg)=nan;
%    u_btrop=squeeze(F);

    [F,n,m,nlr] = read_hycom(fina,finb,'v-vel.');
    F(F>hgg)=0;
    V=F;

%    [F,n,m,l] = read_hycom(fina,finb,'v_btrop');
%    F(F>hgg)=nan;
%    v_btrop=squeeze(F);


    [ZMh,ZZh] = sub_zz_zm(fina,finb,HH);
    dH=abs(diff(ZZh,1));

    nsc=length(SCT);
    for isc=1:nsc
      nin=IPR(isc);
      A = collocate_UTS_section(SCT(isc).I,SCT(isc).J,...
				U,V,T,S,dH,HH,LON,LAT,nin);
      enorm=A.Norm(:,1)+A.Norm(:,2); % need only sign, as norm comp already defined=un
      un=A.Unrm; % normal flow
      Ing=find(enorm<0);
      un(:,Ing)=-un(:,Ing);

      UTS(isc).Hb=A.Hb;
      UTS(isc).Unrm=un;
      UTS(isc).Tnrm=A.Tnrm;
      UTS(isc).Snrm=A.Snrm;
      UTS(isc).dH=A.dH;
      UTS(isc).Segm_dL=A.segmL;
      
      un=UTS(isc).Unrm;
      tn=UTS(isc).Tnrm;
      sn=UTS(isc).Snrm;
      DZ=abs(UTS(isc).dH);
      dx=UTS(isc).Segm_dL;
      
      [DX,dmm]=meshgrid(dx,[1:nlr]);
      
      rhow=sw_dens0(sn,tn);

% Fluxes over sections    
      vf0=nansum(un.*DX.*DZ);
      Vflx=nansum(vf0);
      hf1=nansum(un.*Cp.*rhow.*(tn-Tref1).*DX.*DZ);
      Hflx1=nansum(hf1);
      sf1=nansum(un.*(Sr1-sn)/Sr1.*DX.*DZ);
      Sflx1=nansum(sf1);
      sf2=nansum(un.*(Sr2-sn)/Sr2.*DX.*DZ);
      Sflx2=nansum(sf2);
      
% Update:
      VHFLX(isc).Time=dnmb;
      VHFLX(isc).VolFlx_m3s(cc)=Vflx; % whole section
      VHFLX(isc).HFlx_T1_W(cc)=Hflx1;
      VHFLX(isc).SFlx_S1(cc)=Sflx1;
      VHFLX(isc).SFlx_S2(cc)=Sflx2;
      
% Diangostics:      
      VHFLX(isc).Vfm=VHFLX(isc).Vfm+Vflx;
      VHFLX(isc).Hfm1=VHFLX(isc).Hfm1+Hflx1;
      VHFLX(isc).Sfm1=VHFLX(isc).Sfm1+Sflx1;
      VHFLX(isc).Sfm2=VHFLX(isc).Sfm2+Sflx2;
    end
    

      
    fprintf('\n---------------------------\n');
    for isc=1:nsc
      nm  = VHFLX(isc).Name;
      vf=VHFLX(isc).VolFlx_m3s(cc); % 
      hf1=VHFLX(isc).HFlx_T1_W(cc);
      sf1=VHFLX(isc).SFlx_S1(cc);
      sf2=VHFLX(isc).SFlx_S2(cc);
      
      fprintf('%s Day VolFlx=%5.1fSv, HFlx1=%6.2f TW, SFlx1=%7.3f mSv, SFlx2=%7.3f mSv\n',...
	      nm,vf*1e-6,hf1*1e-12,sf1*1e-3,sf2*1e-3);
    end
    
    for isc=1:nsc
      nm  = VHFLX(isc).Name;
      Vfm = VHFLX(isc).Vfm/cc;
      Hfm1= VHFLX(isc).Hfm1/cc;
      Sfm1= VHFLX(isc).Sfm1/cc;
      Sfm2= VHFLX(isc).Sfm2/cc;
      
      fprintf('%s Mean VolFlx=%5.1fSv, HFlx1=%6.2f TW, SFlx1=%7.3f mSv, SFlx2=%7.3fmSv\n',...
	      nm,Vfm*1e-6,Hfm1*1e-12,Sfm1*1e-3,Sfm2*1e-3);
    end
    
    fprintf('++++++>  Processed 1 record %8.5f min\n\n',toc/60);
      
    if s_mat==1 & mod(cc,60)==0
      fmatout=sprintf('%shycom008_%3.3i_NAtl_fluxes_%4.4i.mat',...
		    pthmat,expt,YR);
      fprintf('Saving %s\n',fmatout);
      save(fmatout,'VHFLX');
    end
  
  end
  
  if s_mat==1
    fmatout=sprintf('%shycom008_%3.3i_NAtl_fluxes_%4.4i.mat',...
		    pthmat,expt,YR);
    fprintf('Saving %s\n',fmatout);
    save(fmatout,'VHFLX');
  end
  
  fprintf(' ========= END OF YEAR %i ======\n\n\n',YR);
  
end
