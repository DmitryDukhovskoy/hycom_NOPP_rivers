% Calculate 
% Tracer fluxes across POP gates
% on Greenland shelf 
%
% The gate coordinates are in the following order: Davis Strait, 
% Cape Farewell, Wide to Narrow, Denmark Strait, Fram Strait, and  Narres Strait. 

addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/Wavelet;
startup;

close all
clear


expt=110;
TV=11;
YR1=2008;
YR2=2008;

s_mat=1;

hgg=1e20;
f_zgrd=0;  % =1 - calculate fluxes from z-grid interpolated U,T,S - less accurate
           % mostly for comparison and validation

pthmat = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_theresa/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';

btx='greenl_Trcr_POPgates.m';

fprintf('arc08-%3.3i Tracer fluxes Greenland Shelf gates %i-%i, save=%i\n',...
	expt,YR1,YR2,s_mat);


ftopo = sprintf('%sdepth_ARCc0.08_%2.2i.nc',pthtopo,TV); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);

GC = sub_greenl_isobath(HH,LON,LAT);
Ig=GC.cntr_Iindx;
Jg=GC.cntr_Jindx;

% Get gates coordinates from POP
flnm=sprintf('%sGates_lat_lons.mat',pthmat);
GT=load(flnm);
pt = [ 0   76    88    97   143   273   287];
NMS{1}='DavisStr';
NMS{2}='Farewell';  % close to OSNAP line
NMS{3}='Tingmiarmiut';  % Wide-to-Narrow SE shelf
NMS{4}='DenmarkStr';
NMS{5}='FramStr';
NMS{6}='NaresStr';
NMS{7}='Fram2';  % control section - straight line
ng=6;

%
% Orientation points to determine 
% +/- fluxes across the sections
% Poisitive is northward or towards the AO
IJPR=[530,748;...
      633,468;...
      709,530;...
      890,623;...
      980,1080;...
      800, 1068;...
      980,1080];
ni=length(IJPR);
for ii=1:ni
  i1=IJPR(ii,1);
  j1=IJPR(ii,2);
  IPR(ii)=sub2ind(size(HH),j1,i1);
end

SCT = struct;
for ii=1:ng
  i1=pt(ii)+1;
  i2=pt(ii+1);
% Start  
  lt1=GT.lat_gts(i1);
  ln1=GT.lon_gts(i1);
  D=distance_spheric_coord(LAT,LON,lt1,ln1);
  [js,is]=find(D==min(min(D)));
% End 
  lt2=GT.lat_gts(i2);
  ln2=GT.lon_gts(i2);
  D=distance_spheric_coord(LAT,LON,lt2,ln2);
  [je,ie]=find(D==min(min(D)));
  
  SCT(ii).Name=NMS{ii};
  SCT(ii).IJ=[is, js; ie, je];

% Correct mismatch with POP grid
  if strncmp(NMS{ii},'Farewell',5)
    SCT(ii).IJ=[612, 421; 634, 421];
  elseif strncmp(NMS{ii},'Tingmiarmiut',5)
    SCT(ii).IJ=[656, 506; 668, 491];
  elseif strncmp(NMS{ii},'NaresStr',5)
    SCT(ii).IJ=[756, 1062; 756, 1073];
  elseif strncmp(NMS{ii},'FramStr',5)
    SCT(ii).IJ=[917, 930; 1073, 950];
  end
  
end
%
% Add straight section for checking
ii=ng+1;
SCT(ii).Name=NMS{ii};
SCT(ii).IJ=[917, 930; 1086, 930];
nsct=length(SCT);

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
%  sub_plot_Greenl_contour(HH,LON,LAT,fn,GC);
  figure(1); clf;
  contour(HH,[0 0],'k');
  hold on;
  contour(HH,[-5000:500:-100],'Color',[0.6 0.6 0.6]);
 
  plot(GC.cntr_Iindx,GC.cntr_Jindx,'b-');
  

  for ip=1:nsct
%    IJ=SCT(ip).IJ;
%    plot([IJ(1,1) IJ(2,1)],[IJ(1,2) IJ(2,2)],...
%	 'Linewidth',2.5,'Color',[1. 0.6 0]);
    IIs=SCT(ip).I;
    JJs=SCT(ip).J;
    plot(IIs,JJs,'-',...
	 'Linewidth',2.5,'Color',[1. 0.3 0]);
    
    Nrm=SCT(ip).Nrm;
    na=length(Nrm);
    for ii=1:na
      x0=0.5*(IIs(ii)+IIs(ii+1));
      y0=0.5*(JJs(ii)+JJs(ii+1));
      x01=x0+2*Nrm(ii,1);
      y01=y0+2*Nrm(ii,2);
      plot([x0 x01],[y0 y01],'c-');
    end
    
    
  end
 
  contour(LON,[-120:10:40],'Color',[0.7 0.7 0.7]);
  contour(LAT,[40:10:88],'Color',[0.7 0.7 0.7]);
%  contour(LAT,[65 65],'Color',[0.8 0.8 0.8]);
  
 
  axis('equal');
  set(gca,'xlim',[450 1100],...
	  'ylim',[350 1200]);
  
  bottom_text(btx,'pwd',1);
  keyboard
end


% Interpolate onto fixed z-grid
% Specify fixed ZZ (interf) and ZM levels:
ZZf = [(0:-1:-10)';(-12:-2:-30)';(-35:-5:-100)';...
       (-110:-10:-1000)';(-1025:-25:-2500)';(-2550:-50:-5000)'];
kzz = length(ZZf);

dZf=abs(diff(ZZf));
ZMf = [];
for ik=1:kzz-1
  ZMf(ik,1)=ZZf(ik)+0.5*(ZZf(ik+1)-ZZf(ik));
end

% Creat dz dx for flux calculation
for isc=1:nsct
  dx=SCT(isc).Segm_dL;
  [DX,DZ]=meshgrid(dx,dZf);
  SCT(isc).DX=DX;
  SCT(isc).DZ=DZ;
end


% Daily Heat, and vol fluxes
VTRCR=SCT;

for YR=YR1:YR2
  yr=YR;
  dE=datenum(yr,12,31);
  dJ1=datenum(yr,1,1);
  ndays=dE-dJ1+1;
  
  cc=0;
  for isc=1:nsct
    VTRCR(isc).Vfm=0;  % track overal means, vol flux
    VTRCR(isc).Trm=0; % tracer flux
  end
  
  for iday=1:7:ndays
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
    nTr=1; % Greenland tracer
    [F,n,m,nlr] = read_hycom(fina,finb,'tracer','r_tracer',nTr);
    F(F>hgg)=nan;
    F(F<0)=0;
    Ctr=F;  % kg/m3

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
      S=Ctr;       % dummy array - nt needed
      A = collocate_UTS_section(SCT(isc).I,SCT(isc).J,...
				U,V,Ctr,S,dH,HH,LON,LAT,nin);
      enorm=A.Norm(:,1)+A.Norm(:,2); % need only sign, as norm comp already defined=un
      un=A.Unrm; % normal flow
      Ing=find(enorm<0);
      if ~isempty(Ing)
        un(:,Ing)=-un(:,Ing);
      end

      UTS(isc).Hb=A.Hb;
      UTS(isc).Unrm=un;
      UTS(isc).Tnrm=A.Tnrm;
      UTS(isc).dH=A.dH;
      UTS(isc).Segm_dL=A.segmL;
      
      un=UTS(isc).Unrm;
      tn=UTS(isc).Tnrm;
      DZ=abs(UTS(isc).dH);
      dx=UTS(isc).Segm_dL;
      
      [DX,dmm]=meshgrid(dx,[1:nlr]);
      
      
      
%      rhow=sw_dens0(sn,tn);

%      Vflx=nansum(nansum(un.*DX.*DZ));
%      Hflx1=nansum(nansum(un.*Cp.*rhow.*(tn-Tref1).*DX.*DZ));
%      Hflx2=nansum(nansum(un.*Cp.*rhow.*(tn-Tref2).*DX.*DZ));
% Fluxes over whole sections    
      vf0=nansum(un.*DX.*DZ);
      Vflx=nansum(vf0);
      TrFlx=nansum(un.*tn.*DX.*DZ);
      tflx=nansum(TrFlx);
      
%      if Vflx*tflx<0
%	fprintf('Vol and Tracer fluxes mismatch sign\n');
%	keyboard
%      end

% Fluxes over Gr Shelf only:
% Need to know where Greenland is
      iGr=SCT(isc).GrCntr_intrcp;
      hb=UTS(isc).Hb;
      if isc~=6
        dh=hb(iGr+5)-hb(iGr-5);
      else
	dh=-1;  % Nares Str, narrow
      end
      
      if dh>0
        Vflx_sh=nansum(vf0(iGr:end));
        TrFlx_sh=nansum(TrFlx(iGr:end));
      else
        Vflx_sh=nansum(vf0(1:iGr));
        TrFlx_sh=nansum(TrFlx(1:iGr));
      end
      
% Update:
      VTRCR(isc).Time=dnmb;
      VTRCR(isc).VolFlx_m3s(cc)=Vflx; % whole section
      VTRCR(isc).TrFlx(cc)=tflx;
      VTRCR(isc).VolFlxGrSh_m3s(cc)=Vflx_sh; % within the contour, shelf
      VTRCR(isc).TrFlx_Sh(cc)=TrFlx_sh; 
      
% Diangostics:      
      VTRCR(isc).Vfm=VTRCR(isc).Vfm+Vflx;
      VTRCR(isc).Trm=VTRCR(isc).Trm+tflx;
    end

    keyboard
    
% Calculate flux from interpolated onto z-grid
% U,T,S fields
% Interpolation by HYCOM layers into Z-grid
    if f_zgrd==1
      UTSZ = sub_sct_interp2z(SCT,U,V,Ctr,S,ZZf,ZZh,HH,dH,DX,DY);    
%
% Calculate Vol and Heat flux
      for isc=1:nsct
	un=UTSZ(isc).NormalU;
	tn=UTSZ(isc).Temp;  % tracer
	dx=SCT(isc).DX;
	dz=abs(SCT(isc).DZ);
	rhow=sw_dens0(sn,tn);
% Should be close to Vflx, Hflx1,2:	
	Vf=nansum(nansum(un.*dx.*dz));
	HTr=nansum(nansum(un.*tn.*dx.*dz));
      end
    
    end        
      
    fprintf('\n---------------------------\n');
    for isc=1:nsc
      nm  = VTRCR(isc).Name;
      vf=VTRCR(isc).VolFlx_m3s(cc); % 
      tr1=VTRCR(isc).TrFlx(cc);
      
      fprintf('%s Day VolFlx=%5.1fSv, TrFlx=%8.2f kg/m3\n',...
	      nm,vf*1e-6,tr1);
    end
    
    for isc=1:nsc
      nm  = VTRCR(isc).Name;
      Vfm = VTRCR(isc).Vfm/cc;
      Trm= VTRCR(isc).Trm/cc;
      
      fprintf('%s Mean VolFlx=%5.1fSv, TrFlxm=%8.2f kg/m3\n',...
	      nm,Vfm*1e-6,Trm);
    end
    
    fprintf('++++++>  Processed 1 record %8.5f min\n\n',toc/60);
      
    if s_mat==1 & mod(cc,20)==0
      fmatout=sprintf('%shycom008_%3.3i_Greenl_Trcr_POPgates_%4.4i.mat',...
		    pthmat,expt,YR);
      fprintf('Saving %s\n',fmatout);
      save(fmatout,'VTRCR');
    end
  
  end
  
  if s_mat==1
    fmatout=sprintf('%shycom008_%3.3i_Greenl_Trcr_POPgates_%4.4i.mat',...
		  pthmat,expt,YR);
    fprintf('Saving %s\n',fmatout);
    save(fmatout,'VTRCR');
  end
  
  fprintf(' ========= END OF YEAR %i ======\n\n\n',YR);
  
end
