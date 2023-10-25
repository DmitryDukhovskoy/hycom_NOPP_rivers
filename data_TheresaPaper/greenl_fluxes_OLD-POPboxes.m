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
YR1=2008;
YR2=2008;

s_mat=1;

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
    SCT(ii).IJ=[656, 506; 668, 491];
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
  
  SCT(ip).Nrm=Nrm;
  SCT(ip).Hbottom=Hbtm;
% Intersection pnt:
% Find the first intercept pont going from the Gr coast 
% in case there are several (when Gr contour and gate segments overlap over several points
%  at the intersection, gowing parallel or perpendicular to the coast)
% This is p-point 
% For closing boxes need UV-contours
  igcst1=SCT(ip).Gr_coast(1);  % check point on the coast
  jgcst1=SCT(ip).Gr_coast(2);
  
  ddc=sqrt((IIs-igcst1).^2+(JJs-jgcst1).^2);

  im=[];
  if ddc(1)>ddc(end)
    for iss=ni:-1:1
      dd=sqrt((Ig-IIs(iss)).^2+(Jg-JJs(iss)).^2);
      mdd=min(dd);
      if mdd==0, break; end
    end
  else
    for iss=1:ni
      dd=sqrt((Ig-IIs(iss)).^2+(Jg-JJs(iss)).^2);
      mdd=min(dd);
      if mdd==0, break; end;
    end
  end
  if abs(mdd)>1e-20
    fprintf('Could not locate GrContour - Gate intercept, section %i: %s\n',ip,SCT(ip).Name);
    keyboard;
  end
  im=iss; 
  jm=find(dd==mdd);

  imn=IIs(im);
  jmn=JJs(im);
  SCT(ip).GrCntr_intrcp=jm;
  SCT(ip).gate_GrIntrcp_I=imn;
  SCT(ip).gate_GrIntrcp_J=jmn;
  SCT(ip).gate_GrIntrcp=im;
  SCT(ip).Segm_dL=dL;
  SCT(ip).Dist_origin=Dst;

  SCT(ip).GrCntr_uvIntrcp=[];  
  SCT(ip).GrCntr_uvI=[];
  SCT(ip).GrCntr_uvJ=[];
end


f_map=1;
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
  for ip=1:nsct
    igr=SCT(ip).gate_GrIntrcp_I;
    jgr=SCT(ip).gate_GrIntrcp_J;
    Icntr=SCT(ip).GrCntr_intrcp;
    icnt=GC.cntr_Iindx(Icntr);
    jcnt=GC.cntr_Jindx(Icntr);

    plot(igr,jgr,'.','Markersize',28,'Color',[0 1 0]);
    plot(icnt,jcnt,'.','Markersize',14,'Color',[1 0 1]);
  end

  axis('equal');
  set(gca,'xlim',[450 1100],...
	  'ylim',[350 1200]);
  
  bottom_text(btx,'pwd',1);

%keyboard
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
BOX(ibx).Pnt=[753, 1068; ...
              539, 677];  % positive direction for each gate section (into the box) 

ibx=ibx+1;
BOX(ibx).Name='DaivsFarw';
BOX(ibx).Sect=[1, 2];    % section #
BOX(ibx).Pnt=[534, 650; ...  % Point inside box to determin + influx 
              615, 415];
ibx=ibx+1;
BOX(ibx).Name='FarwTing';
BOX(ibx).Sect=[2, 3];    % section #
BOX(ibx).Pnt=[619, 423; ... % Point inside box to determin + influx 
              661, 497];
ibx=ibx+1;
BOX(ibx).Name='TingDenm';
BOX(ibx).Sect=[3, 4];    % section #
BOX(ibx).Pnt=[666, 498; ... % Point inside box to determin + influx 
              812, 556];

ibx=ibx+1;
BOX(ibx).Name='DenmFram';
BOX(ibx).Sect=[4, 5];    % section #
BOX(ibx).Pnt=[818, 561;... % Point inside box to determin + influx 
              978, 925];
ibx=ibx+1;
BOX(ibx).Name='FramNares';
BOX(ibx).Sect=[5, 6];    % section #
BOX(ibx).Pnt=[974, 940; ... % Point inside box to determin + influx 
              758, 1075];

Nbx=length(BOX);



for ibx=1:Nbx
% Note that for boxes >1, section 1 = section 2 from previous box
% just flip the orientation of the norms 
  fprintf('Finding corner points to close boxes, ibx=%i %s\n',ibx,BOX(ibx).Name);

  for isct=1:2    
    i1=BOX(ibx).Sect(isct);
    ipbox=BOX(ibx).Pnt(isct,1);
    jpbox=BOX(ibx).Pnt(isct,2);

				igr=SCT(i1).gate_GrIntrcp_I;
				jgr=SCT(i1).gate_GrIntrcp_J;
				Intgt=SCT(i1).gate_GrIntrcp;   % intercept index for the gate section
				Intgr=SCT(i1).GrCntr_intrcp;  % intercept index for the Greenl contour 
				icst=SCT(i1).Gr_coast(1);  % check point on the coast
				jcst=SCT(i1).Gr_coast(2);
				
				IIs=SCT(i1).I; % gate section i pnts
				JJs=SCT(i1).J;

% Section end closest to the coast:
				dd=sqrt((IIs-icst).^2+(JJs-jcst).^2);
				if dd(1)<dd(end)
						st1=1; % Greenl coast at index 1 of the segment
				else
						st1=0;
				end

    if strncmp(SCT(i1).Name,'Nares',5)
      BOX(ibx).intrcp_gate_indx(isct,1)=1;
      BOX(ibx).intrcp_gate_indx(isct,2)=Intgt;
      BOX(ibx).intrcp_grcntr_indx(isct)=1;
      BOX(ibx).Sct(isct).I=IIs(1:Intgt);
      BOX(ibx).Sct(isct).J=JJs(1:Intgt);
      continue
    end 

    if isct==1 & ibx>1
      BOX(ibx).intrcp_gate_indx(isct,1) = BOX(ibx-1).intrcp_gate_indx(2,1);
      BOX(ibx).intrcp_gate_indx(isct,2) = BOX(ibx-1).intrcp_gate_indx(2,2);
      BOX(ibx).Sct(isct).I              = BOX(ibx-1).Sct(2).I;
      BOX(ibx).Sct(isct).J              = BOX(ibx-1).Sct(2).J;
    end

		% Connect corner segments with appropriate U,V points
				[Isct, Igrn] = sub_connect_corner_sgm(IIs,JJs,Ig,Jg,icst,jcst,Intgr,Intgt,ipbox,jpbox,isct,f_map);
    BOX(ibx).intrcp_gate_indx(isct,1)=Isct; % corrected indices section corner point
    if st1==1
      BOX(ibx).intrcp_gate_indx(isct,2)=1;   % section end/start point
      BOX(ibx).Sct(isct).I=IIs(1:Isct);
      BOX(ibx).Sct(isct).J=JJs(1:Isct);
    else
      BOX(ibx).intrcp_gate_indx(isct,2)=length(IIs);
      BOX(ibx).Sct(isct).I=IIs(Isct:end);
      BOX(ibx).Sct(isct).J=JJs(Isct:end);
    end
    BOX(ibx).intrcp_grcntr_indx(isct)=Igrn;
  end
%

end
keyboard


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

% Create dz dx for flux calculation
for isc=1:nsct
  dx=SCT(isc).Segm_dL;
  [DX,DZ]=meshgrid(dx,dZf);
  SCT(isc).DX=DX;
  SCT(isc).DZ=DZ;
end


% Daily Heat, and vol fluxes
VHFLX=SCT;
for isc=1:nsct
  VHFLX(isc).Tref1=Tref1;
  VHFLX(isc).Tref2=Tref2;
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
    VHFLX(isc).Hfm2=0;
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

%    [F,n,m,nlr] = read_hycom(fina,finb,'srfhgt');
%    F(F>hgg)=0;
%    SSH=F/9.806;  % m


    [ZMh,ZZh] = sub_zz_zm(fina,finb,HH);
    dH=abs(diff(ZZh,1));

%
% Colloate T,S -> U,V points
% along the Gr contour
    nin=1; % find norm inside the contour,
    GR = collocate_UTS_section(Ig,Jg,U,V,T,S,dH,HH,LON,LAT,nin);


    nsc=length(SCT);
    for isc=1:nsc
%keyboard
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
      UTS(isc).sgmX=A.sgmX;
      UTS(isc).sgmY=A.sgmY; 
      
      un=UTS(isc).Unrm;
      tn=UTS(isc).Tnrm;
      sn=UTS(isc).Snrm;
      DZ=abs(UTS(isc).dH);
      dx=UTS(isc).Segm_dL;

%
% Find Gr contour - gates intersection
% of UV contours
      dmm=UTS(isc).GrCntr_uvI;
      if isempty(dmm)
        Gx=GR.sgmX;
        Gy=GR.sgmY;
        Sx=A.sgmX;
        Sy=A.sgmY;

      end
      
      [DX,dmm]=meshgrid(dx,[1:nlr]);
      
      rhow=sw_dens0(sn,tn);

%      Vflx=nansum(nansum(un.*DX.*DZ));
%      Hflx1=nansum(nansum(un.*Cp.*rhow.*(tn-Tref1).*DX.*DZ));
%      Hflx2=nansum(nansum(un.*Cp.*rhow.*(tn-Tref2).*DX.*DZ));
% Fluxes over whole sections    
      vf0=nansum(un.*DX.*DZ);
      Vflx=nansum(vf0);
      hf1=nansum(un.*Cp.*rhow.*(tn-Tref1).*DX.*DZ);
      Hflx1=nansum(hf1);
      hf2=nansum(un.*Cp.*rhow.*(tn-Tref2).*DX.*DZ);
      Hflx2=nansum(hf2);

% Fluxes over Gr Shelf only:
% Need to know where Greenland is
      iGr=SCT(isc).GrCntr_intrcp;
      hb=UTS(isc).Hb;
%      dh=hb(iGr+5)-hb(iGr-5);
      if isc~=6
        dh=hb(iGr+5)-hb(iGr-5);
      else
	dh=-1;  % Nares Str, narrow
      end
%keyboard
      if dh>0
	Vflx_sh=nansum(vf0(iGr:end));
	Hflx1_sh=nansum(hf1(iGr:end));
	Hflx2_sh=nansum(hf2(iGr:end));
      else
	Vflx_sh=nansum(vf0(1:iGr));
	Hflx1_sh=nansum(hf1(1:iGr));
	Hflx2_sh=nansum(hf2(1:iGr));
      end
% Update:
      VHFLX(isc).Time=dnmb;
      VHFLX(isc).VolFlx_m3s(cc)=Vflx; % whole section
      VHFLX(isc).HFlx_T1_W(cc)=Hflx1;
      VHFLX(isc).HFlx_T2_W(cc)=Hflx2;
      VHFLX(isc).VolFlxGrSh_m3s(cc)=Vflx_sh; % within the contour, shelf
      VHFLX(isc).HFlxGrSh_T1_W(cc)=Hflx1_sh; 
      VHFLX(isc).HFlxGrSh_T2_W(cc)=Hflx2_sh;
      
% Diangostics:      
      VHFLX(isc).Vfm=VHFLX(isc).Vfm+Vflx;
      VHFLX(isc).Hfm1=VHFLX(isc).Hfm1+Hflx1;
      VHFLX(isc).Hfm2=VHFLX(isc).Hfm2+Hflx2;
    end
    

      
    fprintf('\n---------------------------\n');
    for isc=1:nsc
      nm  = VHFLX(isc).Name;
      vf=VHFLX(isc).VolFlx_m3s(cc); % 
      hf1=VHFLX(isc).HFlx_T1_W(cc);
      hf2=VHFLX(isc).HFlx_T2_W(cc);
      
      fprintf('%s Day VolFlx=%5.1fSv, HFlx1=%6.2f TW, HFlx2=%6.2f TW\n',...
	      nm,vf*1e-6,hf1*1e-12,hf2*1e-12);
    end
    
    for isc=1:nsc
      nm  = VHFLX(isc).Name;
      Vfm = VHFLX(isc).Vfm/cc;
      Hfm1= VHFLX(isc).Hfm1/cc;
      Hfm2= VHFLX(isc).Hfm2/cc;
      
      fprintf('%s Mean VolFlx=%5.1fSv, HFlx1=%6.2f TW, HFlx2=%6.2f TW\n',...
	      nm,Vfm*1e-6,Hfm1*1e-12,Hfm2*1e-12);
    end
    
    fprintf('++++++>  Processed 1 record %8.5f min\n\n',toc/60);
      
    if s_mat==1 & mod(cc,60)==0
      fmatout=sprintf('%shycom008_%3.3i_Greenl_flx_POPgates_%4.4i.mat',...
		    pthmat,expt,YR);
      fprintf('Saving %s\n',fmatout);
      save(fmatout,'VHFLX');
    end
  
  end
  
  if s_mat==1
    fmatout=sprintf('%shycom008_%3.3i_Greenl_flx_POPgates_%4.4i.mat',...
		  pthmat,expt,YR);
    fprintf('Saving %s\n',fmatout);
    save(fmatout,'VHFLX');
  end
  
  fprintf(' ========= END OF YEAR %i ======\n\n\n',YR);
  
end
