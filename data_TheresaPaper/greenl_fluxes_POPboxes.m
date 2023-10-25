% Code uses updated collocation function, provides accurate volume flux calculation
%
% Calculate volume and heat fluxes 
% for the "boxes" on the Greenland shelf
% bounded by the Gr contour and several "gates" 
%
% This code replaces the wrong version greenl_fluxes_POPboxes.m 
%
%
% From Theresa's email 08/13/2019
%
%1) calculate Vol and heat flux across the “gates” - 
% I still need the coordinates of the gates on the shelf
%    The gate lat/lons are: 
%              1. Davis Strait: 67 N
%              2. Cape Farewell: 60 N
%              3. Wide-to-Narrow Gate: 63.5 N
%              4. Denmark Strait: 68.2 N, 30.8W to 65.6 N, 24.2 W 
%              5. Fram Strait: 79.5 N
%              6. Narres Strait (East end): 
%   I think  gates you have defined previously were good and you should only 
% need to define the Wide-to-Narrow Gate. For Gates 1,4,5,6 we need both 
% the full flux across the strait and the flux between the continent and 
% the control volume.  For 2 and 3, only the transport between the continent 
% and contour are needed. I'm attaching a screenshot of the POP/HYCOM 
% gates (POP in blue, HYCOM in red), 
% I've additionally included a .mat file with the lats/lons of the gates in POP.

%2) calculate heat fluxes into and out of boxes on the shelf - monthly
% time series over your whole simulation and daily time series from 2005-2009. 
%3) heat flux across Greenland contour - monthly time series over your 
% whole simulation and daily time series from 2005-2009. 

%I would like to add:
%4) Volume average temperature within the six control volumes - monthly 
% time series over your whole simulation and daily time series from 2005-2009. 
%5) Irminger Basin T/S Samples: I have found the indexes on the HYCOM 
% grid closest to the WOA observations I am using. For each point I just 
% need the vertical profile in the corresponding month 
% and year (also listed in the matrix). 
% HT and KMT are the depth and maximum z-level of the gates respectively. 
%
% I lumped all the gates into one array for convenience, 
% I should have included the indexing vector so you could break it down more easily. 
% pt = [ 0   76    88    97   143   273   287];
% where Davis_St_lat = lat_gt(pt(1)+1:pt(2)) and so on.
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
YR1=2005;
YR2=2005;

s_mat=0;

Cp = 4200; % J/kg K
Tref1= -1.8; % Ref T to calc. H flux
Tref2= 0; % Ref T to calc. H flux
hgg=1e20;
f_zgrd=0;  % =1 - calculate fluxes from z-grid interpolated U,T,S - less accurate
           % mostly for comparison and validation

pthmat = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_theresa/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';

btx='greenl_fluxes_POPboxes.m';

fprintf('arc08-%3.3i Heat and Vol fluxes Greenland Shelf gates %i-%i, save=%i\n',...
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
GC.sgmX=[]; % indices of the contour segments going through UV pnts
GC.sgmY=[]; 


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
NMS{8}='Davis2'; % control section - straight line
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

% Check points on Gr coast for
% defining Gr end of the segments:
% Starting from Daivs Str.
jGR=[583         654
         603         426
         650         513
         776         619
         904         934
         754        1054];

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
    SCT(ii).IJ=[655, 506; 668, 491];
  elseif strncmp(NMS{ii},'NaresStr',5)
    SCT(ii).IJ=[756, 1062; 756, 1083];
  elseif strncmp(NMS{ii},'FramStr',5)
    SCT(ii).IJ=[917, 930; 1073, 950];
  end
  
end
%
% Add straight section for checking
% Fram Str.
ii=ng+1;
SCT(ii).Name=NMS{ii};
SCT(ii).IJ=[917, 930; 1086, 930];
jGR(ii,:)=jGR(5,:);
nsct=length(SCT);

% Davis
ii=nsct+1;
SCT(ii).Name=NMS{ii};
SCT(ii).IJ=[473, 673; 564, 673];
nsct=length(SCT);
jGR(ii,:)=jGR(1,:);


% Change orientation of segments
% segments go from i1 to i2 and i1<i2
% in Nares Str: j1 to j2, j1<j2
P=[0, 1; 1, 0]; % permutation matrix
for ip=1:nsct
  i1=SCT(ip).IJ(1,1);
  i2=SCT(ip).IJ(2,1);
  j1=SCT(ip).IJ(1,2);
  j2=SCT(ip).IJ(2,2);

  if i1>i2
    dmm=SCT(ip).IJ;
    dmm=P*dmm;
    SCT(ip).IJ=dmm;
  elseif i1==i2 & j1>j2
    dmm=SCT(ip).IJ;
    dmm=P*dmm;
    SCT(ip).IJ=dmm;
  end
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
%
% Gr Coast orientation marker
  SCT(ip).Gr_coast(1:2)=jGR(ip,:);
end


%keyboard
% Find intersection with GR contour
%
% Fluxes on the intersection - computed at u/v points
%
%       |         |          |
%       |         |          |
%       |         |          |
%       |         |          |    ^
%  ---------------------------    |
%                                 | Gr Contour going up
%       |         |          |    
%       |  (j,i-1)|   (j,i)  | Intersection point
%       |    *  U --   *     |
%       |         |          |  Flux on contour should be U(j,i)
%       |    |    |    |     |  Flux on gate should stop at V(j,i-1)
%  ----------|---------|------  ---->  Gate going east
%         V(j,i-1)      V(j,i)
%
%
%
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
%keyboard
  SCT(ip).Nrm=Nrm;
  SCT(ip).Hbottom=Hbtm;
  SCT(ip).Segm_dL=dL;
  SCT(ip).Dist_origin=Dst;

  SCT(ip).GrCntr_uvIntrcp=[];  
  SCT(ip).GrCntr_uvI=[];
  SCT(ip).GrCntr_uvJ=[];
end


f_map=0;
if f_map==1
  fprintf('Drawing map with segments\n');
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
  end
  
% Plot coast markers
  for ip=1:nsct
    i1=SCT(ip).Gr_coast(1);
    j1=SCT(ip).Gr_coast(2);

    plot(i1,j1,'.','Markersize',18,'Color',[0 0.5 0.8]);
    text(i1,j1+ip,SCT(ip).Name);
  end

% Plot intersection points:
%  for ip=1:nsct
%    igr=SCT(ip).gate_GrIntrcp_I;
%    jgr=SCT(ip).gate_GrIntrcp_J;
%    Icntr=SCT(ip).GrCntr_intrcp;
%    icnt=GC.cntr_Iindx(Icntr);
%    jcnt=GC.cntr_Jindx(Icntr);
%
%    plot(igr,jgr,'.','Markersize',28,'Color',[0 1 0]);
%    plot(icnt,jcnt,'.','Markersize',14,'Color',[1 0 1]);
%  end

  axis('equal');
  set(gca,'xlim',[450 1100],...
	  'ylim',[350 1200]);
  
  bottom_text(btx,'pwd',1);

keyboard
end

CLR=[0 0.4 0.8; ...
     0.8 0.4 0; ...
     1 0.2 0; ...
     0 1 0; ...
     0.8 0 0.6; ...
     0.8 1 0];


% Construct boxes
ibx=1;
BOX(ibx).Name='NaresDavis';
BOX(ibx).Sect=[ 6, 1];    % sections # 
BOX(ibx).Pnt=[727, 1063; ...
              583, 800];  % positive direction for each gate section (into the box) 

ibx=ibx+1;
BOX(ibx).Name='DaivsFarw';
BOX(ibx).Sect=[1, 2];    % section #
BOX(ibx).Pnt=[529, 574; ...  % Point inside box to determin + influx 
              609, 404];
ibx=ibx+1;
BOX(ibx).Name='FarwTing';
BOX(ibx).Sect=[2, 3];    % section #
BOX(ibx).Pnt=[637, 498; ... % Point inside box to determin + influx 
              614, 440];
ibx=ibx+1;
BOX(ibx).Name='TingDenm';
BOX(ibx).Sect=[3, 4];    % section #
BOX(ibx).Pnt=[786, 581; ... % Point inside box to determin + influx 
              669, 507];

ibx=ibx+1;
BOX(ibx).Name='DenmFram';
BOX(ibx).Sect=[4, 5];    % section #
BOX(ibx).Pnt=[988, 847;... % Point inside box to determin + influx 
              896, 653];
ibx=ibx+1;
BOX(ibx).Name='FramNares';
BOX(ibx).Sect=[5, 6];    % section #
BOX(ibx).Pnt=[924, 1030; ... % Point inside box to determin + influx 
              882, 1066];

nbx=length(BOX);


% Find contour through UV pnts for Greenland contour
nin=1; % positive inside the contour
UVGR = sub_UVpnts_contour(Ig,Jg,nin,HH);
dL=sub_segm_dL(UVGR.sgmX,UVGR.sgmY,LON,LAT);
GC.Norm=UVGR.Norm;
GC.sgmX=UVGR.sgmX;
GC.sgmY=UVGR.sgmY;
GC.Ictr=UVGR.Ictr;
GC.IJ_indx=UVGR.gridIndx_IJ; % grid indices corresponding to UV segments
GC.IJ_adj_indx=UVGR.adjIndx_I1J1; % adjacent grid point index - for collocating dH to UV pnts
GC.segm_dL=dL;

for ibx=1:nbx
% Note that for boxes >1, section 1 = section 2 from previous box
% just flip the orientation of the norms 
  fprintf('Finding UV-point segments, box=%i %s\n',ibx,BOX(ibx).Name);

  for isct=1:2    
    i1=BOX(ibx).Sect(isct);
    ipbox=BOX(ibx).Pnt(isct,1);
    jpbox=BOX(ibx).Pnt(isct,2);

    BOX(ibx).S(isct).Name=SCT(i1).Name;

				IIs=SCT(i1).I; % gate section i pnts
				JJs=SCT(i1).J;

    if isct==1 & ibx>1
						BOX(ibx).S(isct).Norm=-BOX(ibx-1).S(2).Norm;
						BOX(ibx).S(isct).sgmX=BOX(ibx-1).S(2).sgmX;
						BOX(ibx).S(isct).sgmY=BOX(ibx-1).S(2).sgmY;
      BOX(ibx).S(isct).IJ_indx=BOX(ibx-1).S(2).IJ_indx; % grid indices corresponding to UV segments
      BOX(ibx).S(isct).adjIndx_I1J1=BOX(ibx-1).S(2).adjIndx_I1J1; % adjacent grid point index - for collocating dH to UV pnts
      BOX(ibx).S(isct).segm_dL=BOX(ibx-1).S(2).segm_dL; 
      continue
    end
    if isct==2 & ibx==nbx
						BOX(ibx).S(isct).Norm=-BOX(1).S(1).Norm;
						BOX(ibx).S(isct).sgmX=BOX(1).S(1).sgmX;
						BOX(ibx).S(isct).sgmY=BOX(1).S(1).sgmY;
      BOX(ibx).S(isct).IJ_indx=BOX(1).S(1).IJ_indx; % grid indices corresponding to UV segments
      BOX(ibx).S(isct).adjIndx_I1J1=BOX(1).S(1).adjIndx_I1J1; % adjacent grid point index - for collocating dH to UV pnts
      BOX(ibx).S(isct).segm_dL=BOX(1).S(1).segm_dL;
      continue
    end
%
% Find segments of the UV contour for gate section
    nin=sub2ind(size(HH),jpbox,ipbox);
    UVGR = sub_UVpnts_contour(IIs,JJs,nin,HH);
    dL=sub_segm_dL(UVGR.sgmX,UVGR.sgmY,LON,LAT);

    BOX(ibx).S(isct).Norm=UVGR.Norm;
    BOX(ibx).S(isct).sgmX=UVGR.sgmX;
    BOX(ibx).S(isct).sgmY=UVGR.sgmY;
    BOX(ibx).S(isct).IJ_indx=UVGR.gridIndx_IJ; % grid indices corresponding to UV segments
    BOX(ibx).S(isct).adjIndx_I1J1=UVGR.adjIndx_I1J1; % adjacent grid point index - for collocating dH to UV pnts
    BOX(ibx).S(isct).segm_dL=dL;

  end
%

end
%keyboard

% Find intersection points and close boxes
for ibx=1:nbx
  for isct=1:2;
				i1=BOX(ibx).Sect(isct);

    fprintf('Finding intersection, Box %i sect=%s\n',ibx,SCT(i1).Name);

				ipbox=BOX(ibx).Pnt(isct,1); % point in the box
				jpbox=BOX(ibx).Pnt(isct,2);

    IGate=BOX(ibx).S(isct).IJ_indx(:,1);  % gate section points
    JGate=BOX(ibx).S(isct).IJ_indx(:,2);
	%			igr=SCT(i1).gate_GrIntrcp_I;
	%			jgr=SCT(i1).gate_GrIntrcp_J;
	%			IGt=SCT(i1).gate_GrIntrcp;   % intercept index for the gate section
	%			IGr=SCT(i1).GrCntr_intrcp;  % intercept index for the Greenl contour 
				icst=SCT(i1).Gr_coast(1);  % check point on the coast
				jcst=SCT(i1).Gr_coast(2);
				IIs=BOX(ibx).S(isct).sgmX;
				JJs=BOX(ibx).S(isct).sgmY;
				Igs=GC.sgmX;
				Jgs=GC.sgmY;
    IGrnl=GC.IJ_indx(:,1);
    JGrnl=GC.IJ_indx(:,2);
			
				[IxGr,IxGt,IxGt0] = sub_intrcpt(IIs,JJs,IGate,JGate,Igs,Jgs,IGrnl,JGrnl,icst,jcst,ipbox,jpbox);
    BOX(ibx).S(isct).Intrcp_Greenl=IxGr;
    BOX(ibx).S(isct).Intrcp_Gate=IxGt;
    BOX(ibx).S(isct).Start_Gate=IxGt0;
  end
end

f_box=0;
if f_box==1
% Plot all boxes and sections and norms
% for each box
  fprintf('Drawing map with boxes\n');
  nf=15;
  sub_check_boxes(nf,HH,GC,BOX);
  keyboard
end


for YR=YR1:YR2
  yr=YR;
  dE=datenum(yr,12,31);
  dJ1=datenum(yr,1,1);
  ndays=dE-dJ1+1;
  
  cc=0;
  
% Daily Heat, and vol fluxes
		VHFLX=struct;
		for ibx=1:nbx
				for isct=1:2
						VHFLX(ibx).S(isct).Tref1=Tref1;
						VHFLX(ibx).S(isct).Tref2=Tref2;
						VHFLX(ibx).S(isct).Vfm=0;
						VHFLX(ibx).S(isct).Hf1m=0;
						VHFLX(ibx).S(isct).Hf2m=0;
						VHFLX(ibx).S(isct).VfmGr=0;
						VHFLX(ibx).S(isct).Hf1mGr=0;
						VHFLX(ibx).S(isct).Hf2mGr=0;
				end
		end
		GVHFLX=struct;
		GVHFLX.Tref1 = Tref1;
		GVHFLX.Tref2 = Tref2;
		GVHFLX.Vfm   = 0;
		GVHFLX.Hf1m  = 0; 
		GVHFLX.Hf2m  = 0;

  for iday=1:ndays
%  for iday=282:282
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

%    [F,n,m,nlr] = read_hycom(fina,finb,'srfhgt');
%    F(F>hgg)=0;
%    SSH=F/9.806;  % m


    [ZMh,ZZh] = sub_zz_zm(fina,finb,HH);
    dH=abs(diff(ZZh,1));

%  =========================================
% Calculate fluxes through Gr Contour
%  =========================================
    chgr=0;
    if chgr==1
      figure(17); clf; hold on;
      contour(HH,[0 0],'k');
      axis('equal');
      set(gca,'xlim',[420 1000],...
              'ylim',[350 1100]);
    end

    iGRs=GC.sgmX;
    jGRs=GC.sgmY;
    nns=length(iGRs(:,1));
    fprintf('Fluxes through Gr Contour ...\n');
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
% Collocate
						if y1==y2  % horiz segment, V flux
								CLC=sub_collocate2uv(V,T,S,dH,HH,i0,j0,i1,j1);
						else
								CLC=sub_collocate2uv(U,T,S,dH,HH,i0,j0,i1,j1);
						end
						In=find(CLC.dHn<1e-3);
						CLC.Tn(In)=nan;
						CLC.Sn(In)=nan;
						CLC.Un(In)=nan;

      snrm=sign(xnrm+ynrm); % one is 0,  need only direction
      GUTS.Hb(j)=-CLC.Hn;
      GUTS.Unrm(:,j)=snrm*CLC.Un;
      GUTS.Tnrm(:,j)=CLC.Tn;
      GUTS.Snrm(:,j)=CLC.Sn;
      GUTS.dH(:,j)=CLC.dHn;

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
    tn=GUTS.Tnrm;
    sn=GUTS.Snrm;
    DZ=abs(GUTS.dH);
    dx=GC.segm_dL;
    [DX,dmm]=meshgrid(dx,[1:nlr]);
    rhow=sw_dens0(sn,tn);

% Fluxes over whole sections    
    vfG0=nansum(un.*DX.*DZ);
    hfG1=nansum(un.*Cp.*rhow.*(tn-Tref1).*DX.*DZ);
    hfG2=nansum(un.*Cp.*rhow.*(tn-Tref2).*DX.*DZ);
    GVHFLX.Time(cc)=dnmb;
    GVHFLX.Vflx(cc)=nansum(vfG0); % overall vol flux
    GVHFLX.Hflx1(cc)=nansum(hfG1);  % overall heat flux
    GVHFLX.Hflx2(cc)=nansum(hfG2);  % overall heat flux ref T2
%
%   Fluxes pointwise across the Gr. contour:
    GVHFLX.VflxPnts(cc,:)=vfG0;
    GVHFLX.Hflx1Pnts(cc,:)=hfG1;
    GVHFLX.Hflx2Pnts(cc,:)=hfG2;
% Mean:
    GVHFLX.Vfm = GVHFLX.Vfm+nansum(vfG0);
    GVHFLX.Hf1m = GVHFLX.Hf1m+nansum(hfG1);
    GVHFLX.Hf2m = GVHFLX.Hf2m+nansum(hfG2);

%  =============================================================
% Calculate fluxes through sections/ Gr contour by boxes

   for ibx=1:nbx
     UTS=struct;
     for isct=1  % 2nd section = -1* 1st section of next box
							iGTs=BOX(ibx).S(isct).sgmX;   % i indices of the uv segments start/end points
							jGTs=BOX(ibx).S(isct).sgmY;   % j indices of the uv segments
							Igtx=BOX(ibx).S(isct).Intrcp_Gate;  % end/intercept index of gate sections
							Igt0=BOX(ibx).S(isct).Start_Gate;   % start/coast index of gate section
							IJp=BOX(ibx).S(isct).IJ_indx;               % grid index corrsponding to the segment
							IJadj=BOX(ibx).S(isct).adjIndx_I1J1;        % corresponding adjacent grid cell for collocation
							Nrm=BOX(ibx).S(isct).Norm;                  % unit norm in the box
       segm_dL=BOX(ibx).S(isct).segm_dL;           % segment length, m

							if Igt0>Igtx                
									ii1=Igtx;                                 
									ii2=Igt0;                                 
							else                        
									ii1=Igt0;                                 
									ii2=Igtx;                                 
							end                         
 
       nns=length(iGTs(:,1));  
							for j=1:nns
									x1=iGTs(j,1);
									x2=iGTs(j,2);
									y1=jGTs(j,1);
									y2=jGTs(j,2);
         xnrm=Nrm(j,1);
         ynrm=Nrm(j,2);
%
% Grid cell and adjacent cell
         i0=IJp(j,1);
         j0=IJp(j,2);
         i1=IJadj(j,1);
         j1=IJadj(j,2);
% Collocate
         if y1==y2  % horiz segment, V flux
           CLC=sub_collocate2uv(V,T,S,dH,HH,i0,j0,i1,j1);
         else
           CLC=sub_collocate2uv(U,T,S,dH,HH,i0,j0,i1,j1);
         end
         In=find(CLC.dHn<1e-3);
         CLC.Tn(In)=nan;
         CLC.Sn(In)=nan;
         CLC.Un(In)=nan;

         snrm=sign(xnrm+ynrm); % one is 0,  need only direction
									UTS.Hb(j)=-CLC.Hn;
									UTS.Unrm(:,j)=snrm*CLC.Un;
									UTS.Tnrm(:,j)=CLC.Tn;
									UTS.Snrm(:,j)=CLC.Sn;
									UTS.dH(:,j)=CLC.dHn;

         fchck=0;
         if fchck==1
%											figure(16); clf; hold on
											plot([x1 x2],[y1 y2],'-');
											plot(i0,j0,'r*');
											plot(i1,j1,'g*');
											xp=0.5*(x1+x2);
											yp=0.5*(y1+y2);
											plot([xp xp+xnrm],[yp yp+ynrm],'g-');
        end

      end  % segment in 1 section
      un=UTS.Unrm;
      tn=UTS.Tnrm;
      sn=UTS.Snrm;
      DZ=abs(UTS.dH);
      dx=segm_dL;

      [DX,dmm]=meshgrid(dx,[1:nlr]);
      
      rhow=sw_dens0(sn,tn);

% Fluxes over whole sections    
      vf0=nansum(un.*DX.*DZ);
      Vflx=nansum(vf0);
      hf1=nansum(un.*Cp.*rhow.*(tn-Tref1).*DX.*DZ);
      Hflx1=nansum(hf1);
      hf2=nansum(un.*Cp.*rhow.*(tn-Tref2).*DX.*DZ);
      Hflx2=nansum(hf2);

% Fluxes over Gr Shelf only, boxes
     	Vflx_sh=nansum(vf0(ii1:ii2));
						Hflx1_sh=nansum(hf1(ii1:ii2));
						Hflx2_sh=nansum(hf2(ii1:ii2));
%
% Greenland contour: 
% by default, fluxes through Gr. contour are + into the box (on Gr Shelf)
					 igx1=BOX(ibx).S(1).Intrcp_Greenl;  % start/end indices of Greenl contour
					 igx2=BOX(ibx).S(2).Intrcp_Greenl;  % start/end indices of Greenl contour
      if ibx==1, % there is discont in indices in Nares Str, contour cut line
        igx1=1; % start from 1 - no flux over the land
      end
  
      if igx1>igx2
        dmm=igx1;
        igx1=igx2;
        igx2=dmm;
      end

      VHFLX(ibx).S(isct).Vflx_Gr(cc)=nansum(vfG0(igx1:igx2));
      VHFLX(ibx).S(isct).Hflx1_Gr(cc)=nansum(hfG1(igx1:igx2));
      VHFLX(ibx).S(isct).Hflx2_Gr(cc)=nansum(hfG2(igx1:igx2));

% Check
%      GrVflx = vfG0(igx1:igx2);
%      GrHf1  = hfG1(igx1:igx2);
%      GrHf2  = hfG2(igx1:igx2);
%
%      VHf = GrVflx.*GrHf1;
%      I = find(VHf<0);
%
%      if ~isempty(I) & ibx==3; 
%        fprintf('Check Vol flux and Hflux (-1.8) - different signs\n');
%        keyboard
%      end

% Update:
      VHFLX(ibx).Time(cc)=dnmb;
      VHFLX(ibx).S(isct).VolFlx_m3s(cc)=Vflx; % whole section
      VHFLX(ibx).S(isct).HFlx_T1_W(cc)=Hflx1;
      VHFLX(ibx).S(isct).HFlx_T2_W(cc)=Hflx2;
      VHFLX(ibx).S(isct).VolFlxGrSh_m3s(cc)=Vflx_sh; % within the contour, shelf
      VHFLX(ibx).S(isct).HFlxGrSh_T1_W(cc)=Hflx1_sh; 
      VHFLX(ibx).S(isct).HFlxGrSh_T2_W(cc)=Hflx2_sh;
      
% Diagnostics:      
      VHFLX(ibx).S(isct).Vfm=VHFLX(ibx).S(isct).Vfm+Vflx_sh;
      VHFLX(ibx).S(isct).Hf1m=VHFLX(ibx).S(isct).Hf1m+Hflx1_sh;
      VHFLX(ibx).S(isct).Hf2m=VHFLX(ibx).S(isct).Hf2m+Hflx2_sh;
      VHFLX(ibx).S(isct).VfmGr=VHFLX(ibx).S(isct).VfmGr+nansum(vfG0(igx1:igx2));
      VHFLX(ibx).S(isct).Hf1mGr=VHFLX(ibx).S(isct).Hf1mGr+nansum(hfG1(igx1:igx2));
      VHFLX(ibx).S(isct).Hf2mGr=VHFLX(ibx).S(isct).Hf2mGr+nansum(hfG2(igx1:igx2));
    end   % sections 
    
  end % boxes
    
%keyboard
% Diagnostics
  
    fprintf('\n---------------------------\n');
% Greenl Contour fluxes: total
%    VfGall   = GVHFLX.Vflx;
%    Hf1Gall  = GVHFLX.Hflx1;
%    Hf2Gall  = GVHFLX.Hflx2;
    VfmG0  = GVHFLX.Vfm/cc;
    Hf1mG0 = GVHFLX.Hf1m/cc;
    Hf2mG0 = GVHFLX.Hf2m/cc;

    for ibx=1:nbx
      
      nm    = BOX(ibx).Name;
      isc1  = BOX(ibx).Sect(1);
      isc2  = BOX(ibx).Sect(2);
      sctn1 = SCT(isc1).Name;
      sctn2 = SCT(isc2).Name;
% Section1:
      vfs1  = VHFLX(ibx).S(1).VolFlxGrSh_m3s(cc); % 
      hf1s1 = VHFLX(ibx).S(1).HFlxGrSh_T1_W(cc);
      hf2s1 = VHFLX(ibx).S(1).HFlxGrSh_T2_W(cc);
% Mean:      
      Vfm1  = VHFLX(ibx).S(1).Vfm/cc;  % mean vol flux through section 1
      Hf1m1 = VHFLX(ibx).S(1).Hf1m/cc; % heat flx T1
      Hf2m1 = VHFLX(ibx).S(1).Hf2m/cc; % heat flx T2
      
      ibx0=ibx+1;
      if (ibx==6) ibx0=1; end
% Section2:
      vfs2  = -VHFLX(ibx0).S(1).VolFlxGrSh_m3s(cc); % 
      hf1s2 = -VHFLX(ibx0).S(1).HFlxGrSh_T1_W(cc);
      hf2s2 = -VHFLX(ibx0).S(1).HFlxGrSh_T2_W(cc);
% Mean:      
      Vfm2  = -VHFLX(ibx0).S(1).Vfm/cc;  % mean vol flux through section 1
      Hf1m2 = -VHFLX(ibx0).S(1).Hf1m/cc; % heat flx T1
      Hf2m2 = -VHFLX(ibx0).S(1).Hf2m/cc; % heat flx T2

% Greenl contour flux in the box:
      VfG  = VHFLX(ibx).S(1).Vflx_Gr(cc); 
      hf1G = VHFLX(ibx).S(1).Hflx1_Gr(cc);
      hf2G = VHFLX(ibx).S(1).Hflx2_Gr(cc);

% Gr flux mean:
      VfGm = VHFLX(ibx).S(1).VfmGr/cc;
      Hf1Gm= VHFLX(ibx).S(1).Hf1mGr/cc;
      Hf2Gm= VHFLX(ibx).S(1).Hf2mGr/cc;

% Net fluxes: daily
      VFlx=vfs1+vfs2+VfG;
      HFlx1=hf1s1+hf1s2+hf1G; 
      HFlx2=hf2s1+hf2s2+hf2G; 
% Mean
      VflxM=(Vfm1+Vfm2+VfGm);
      HFlx1M=(Hf1m1+Hf1m2+Hf1Gm);
      HFlx2M=(Hf2m1+Hf2m2+Hf2Gm);
% Error relative (assuming mean net flux=0):
      Verr=VflxM/(abs(Vfm1)+abs(Vfm2)+abs(VfGm)); 

      fprintf('%s Daily fluxes (+ into box): Section 1, %s:\n',nm, sctn1);
      fprintf('Vol=%5.1fSv, Heat1=%6.2f TW, Heat2=%6.2f TW\n',...
	      vfs1*1e-6,hf1s1*1e-12,hf2s1*1e-12);
      fprintf('Section 2, %s:\n',sctn2);
      fprintf('Vol=%5.1fSv, Heat1=%6.2f TW, Heat2=%6.2f TW\n',...
	      vfs2*1e-6,hf1s2*1e-12,hf2s2*1e-12);
      fprintf('Greenl Contour into the box:\n')
      fprintf('Vol=%5.1fSv, Heat1=%6.2f TW, Heat2=%6.2f TW\n',...
	      VfG*1e-6,hf1G*1e-12,hf2G*1e-12);
      fprintf('==== Net fluxes into the box:\n');
      fprintf('Vol=%5.1fSv, Heat1=%6.2f TW, Heat2=%6.2f TW\n',...
	      VFlx*1e-6,HFlx1*1e-12,HFlx2*1e-12);

      fprintf('\n>>>>>>> Mean fluxes, # records=%s \n',cc);
      fprintf('    Section 1, %s: \n',sctn1);
      fprintf('Vol=%5.1fSv, Heat1=%6.2f TW, Heat2=%6.2f TW\n',...
	     Vfm1*1e-6,Hf1m1*1e-12,Hf2m1*1e-12);
      fprintf('    Section 2, %s:\n', sctn2);
      fprintf('Vol=%5.1fSv, Heat1=%6.2f TW, Heat2=%6.2f TW\n',...
	     Vfm2*1e-6,Hf1m2*1e-12,Hf2m2*1e-12);
      fprintf('  Greenl Contour:\n');
      fprintf('Vol=%5.1fSv, Heat1=%6.2f TW, Heat2=%6.2f TW\n',...
	     VfGm*1e-6,Hf1Gm*1e-12,Hf2Gm*1e-12);

      fprintf('======  Net Mean Fluxes #records=%i\n',cc);
      fprintf('Vol=%5.1fSv, Heat1=%6.2f TW, Heat2=%6.2f TW\n',...
	     VflxM*1e-6,HFlx1M*1e-12,HFlx2M*1e-12);
      fprintf('Relative error: %8.5f\n',Verr);
      fprintf(' ............   \n');
      fprintf('Overall GrContour flux and rel.error: %8.4f Sv %8.4f%%\n',...
              nansum(vfG0)*1e-6, nansum(vfG0)/nansum(abs(vfG0))*100);
      fprintf('Overall heat flux T1= %4.1f, T2=%4.1f: %8.6fTW, %8.6fTW\n',...
              Tref1, Tref2, GVHFLX.Hflx1(cc)*1e-12, GVHFLX.Hflx2(cc)*1e-12);
      fprintf('  ==================== \n\n');

    end
    
    
    fprintf('++++++>  Processed 1 record %8.5f min\n\n',toc/60);

%keyboard
      
    if s_mat==1 & mod(cc,60)==0
      fmatout=sprintf('%shycom008_%3.3i_Greenl_flx_POPboxes_%4.4i.mat',...
		    pthmat,expt,YR);
      fprintf('Saving %s\n',fmatout);
      save(fmatout,'VHFLX','BOX','GVHFLX');
    end
  
  end
end  

if s_mat==1
		fmatout=sprintf('%shycom008_%3.3i_Greenl_flx_POPboxes_%4.4i.mat',...
		pthmat,expt,YR);
		fprintf('Saving %s\n',fmatout);
		save(fmatout,'VHFLX','BOX','GVHFLX');
end

fprintf(' ========= END OF YEAR %i ======\n\n\n',YR);


