% This code is similar to ocn_hflx_greenl008.m
% but calculates daily (instead of monthly) mean depth-integrated 
% values 
%
% Calculate ocean vol/heat flux to Greenland
% across specified contour -
% isobath around Greenland
%
% Also calculates volume-mean T within the boxes
% on the shelf bounded by Gr contour and gates
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/seawater;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

dday = 1;
uadj  = 1; % adjust vol flux to be balanced along the contour

YR1 = 2008;
YR2 = 2008;


Cp = 4200; % J/kg K
%Tref= -273.15; % Ref T to calc. H flux
Tref1= -1.8; % Ref T to calc. H flux
Tref2=0;

expt = 110;

s_mat = 1; % =0 - extract data, no save; =1 extract & save, =2 - load saved


rg=9806;  % convert pressure to depth, m
hgg=1e20; 

plr=0; % highlight this interface
btx = 'vhflx_greenl_cntr008.m';

fprintf('Oceanic Heat Flux, Greenland Section, %i-%i\n',YR1,YR2);

regn = 'ARCc0.08';
expt = 110; % experiment without runoff
%expt = 112;  % epxeriment with Greenland runoff and monthly Arctic rivers
pthfig=sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/%s/%3.3i/fig_green_xsct/',...
		  regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
%pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
pthmat  = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_theresa/';

ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DXh,DYh]=sub_dx_dy(LON,LAT);
Acell=DXh.*DYh;

GC = sub_greenl_isobath(HH,LON,LAT);
Hs = GC.Hbottom; % bottom along section
nh = length(Hs);
dL = diff(GC.Distance_m);
dL(nh)=dL(nh-1); % over land anyway, dL is unimportant
GC.dL = dL;

% ================================================
% Plot Greenland map and the contour
% ================================================
f_pltgr=0;
if f_pltgr==1
  fn=10;
  sub_plot_Greenl_contour(HH,LON,LAT,fn,GC);
  
  bottom_text(btx,'pwd',1);
end

% 
if s_mat==1
  fprintf('Mat file will be saved %s\n',pthmat);
end



% ================================================
HFLX = struct;
HFLX.Info = 'Heat & Vol. Fluxes across 800m isobath around Greenland';
HFLX.Tref1     = Tref1; 
HFLX.Tref2     = Tref2; 
HFLX.GrCntr_II = GC.cntr_Iindx;
HFLX.GrCntr_JJ = GC.cntr_Jindx;
HFLX.Hbottom   = GC.Hbottom(:);
HFLX.DistCntr  = GC.Distance_m(:);

% Find intersection points with POP gates on Gr Shelf:
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
ng=6;

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

% Find interception points
% of the Gr contour and Greenland shelf gates
II=HFLX.GrCntr_II;
JJ=HFLX.GrCntr_JJ;
for isc=1:nsct
  Is=SCT(isc).I;
  Js=SCT(isc).J;
  nsgm=length(Is);
  clear D
  for isg=1:nsgm
    is0=Is(isg);
    js0=Js(isg);
    d=sqrt((II-is0).^2+(JJ-js0).^2);
    icr=find(d==min(d),1);
    D(isg,1)=min(d);
    D(isg,2)=icr;  % contour index
  end
  imn=find(D==min(D),1);  % gate index
  HFLX.GrCntrIntercept(isc,1)=D(imn,1);
  HFLX.GrCntrIntercept(isc,2)=D(imn,2);
  HFLX.GateIntercept(isc,1)=imn;
end

% Find points inside the shelf boxes
[IHX,JHX]=meshgrid([1:nn],[1:mm]);
for isc=1:nsct
  fprintf('searching in points, %i\n',isc);
  isc1=isc;
  if isc1==nsct;
    isc2=1;
  else
    isc2=isc+1;
  end

  Is1=SCT(isc1).I;
  Js1=SCT(isc1).J;
  Is2=SCT(isc2).I;
  Js2=SCT(isc2).J;
  Is1=Is1(:);
  Is2=Is2(:);
  Js1=Js1(:);
  Js2=Js2(:);
  icp1=HFLX.GateIntercept(isc1);
  icp2=HFLX.GateIntercept(isc2);
  icr1=HFLX.GrCntrIntercept(isc1,2);
  icr2=HFLX.GrCntrIntercept(isc2,2);
% Check if start or end point is inside the contour
% i.e. if the gate starts from Greenland or 
% ends on Greenland
  inp1s=inpolygon(Is1(1),Js1(1),II,JJ);
  inp1e=inpolygon(Is1(end),Js1(end),II,JJ);
  inp2s=inpolygon(Is2(1),Js2(1),II,JJ);
  inp2e=inpolygon(Is2(end),Js2(end),II,JJ);
  
  if inp1s
    Igate1=Is1(1:icp1);
    Jgate1=Js1(1:icp1);
  else
    Igate1=Is1(end:-1:icp1);
    Jgate1=Js1(end:-1:icp1);
  end
  
  if inp2s
    Igate2=Is2(icp2:-1:1);
    Jgate2=Js2(icp2:-1:1);
  else
    Igate2=Is2(icp2:end);
    Jgate2=Js2(icp2:end);
  end
  
  ic0=606;
  jc0=687;
  if icr1>icr2
    Ibx=[Igate1; II(icr1:end); II(1:icr2); Igate2; ic0];
    Jbx=[Jgate1; JJ(icr1:end); JJ(1:icr2); Jgate2; jc0];
  else
    Ibx=[Igate1; II(icr1:icr2); Igate2; ic0];
    Jbx=[Jgate1; JJ(icr1:icr2); Jgate2; jc0];
  end
  
  INP=inpolygon(IHX,JHX,Ibx,Jbx);
  IN=find(INP==1 & HH<0);
  SCT(isc).Inpoints=IN;
  SCT(isc).Ibx=Ibx;
  SCT(isc).Jbx=Jbx;
  
end

f_chck=0;
if f_chck==1
  figure(10); clf;
  contour(HH,[0 0],'k');
  hold on
  axis('equal');
  set(gca,'xlim',[440 1080],...
	  'ylim',[350 1100]);
  plot(II,JJ,'r.-');
  for isc=1:nsct
    IN=SCT(isc).Inpoints;
    Ibx=SCT(isc).Ibx;
    Jbx=SCT(isc).Jbx;
    plot(IHX(IN),JHX(IN),'.','Color',[0 1-isc/nsct isc/nsct]);
    plot(Ibx,Jbx,'-');
  end
  
  title('POP boxes and gates on Greenland Shelf');
  bottom_text(btx,'pwd',1);
  
end
%keyboard

% save contour information
% with intercept points
fgrd=sprintf('%shycom008_GrContourFlux_info.mat',pthmat);
fprintf('Saving %s\n',fgrd);
GRCNTR=HFLX;
save(fgrd,'GRCNTR');



mold  = 0;
for iyr=YR1:YR2
  cc=0;
  TM = [];
  yr = iyr;
  Vmn=0;
  Hmn1=0;
  Hmn2=0;

  SCT(isc).Tbox=[];
  SCT(isc).Sbox=[];
  SCT(isc).TM=[];

  HFLX=struct;
  fmat = sprintf('%s%3.3i_GreenlCntr_HVflx_daily_%i.mat',...
		   pthmat,expt,iyr);
  
  for iday = 1:dday:366
    pthbin = sprintf('/nexsan/archive/ARCc0.08_%3.3i/data/%4.4i/',expt,yr);
    fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
    finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);
  
    cc   = cc+1;
    j1d  = datenum(yr,1,1);
    dnmb = j1d+iday-1;
    DV   = datevec(dnmb);
    imo  = DV(2);
    TM(cc,1) = dnmb;

    
    fprintf('Reading %4.4i/%2.2i/%2.2i: %s\n',DV(1:3),fina);
    if ~exist(fina,'file');
      fprintf('Not found %s\n',fina);
      continue;
    end
    
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

    fld='thknss';
    [F,n,m,l] = read_hycom(fina,finb,fld);
%    [F,n,m,l] = read_hycom(fina,finb,fld,'r_layer',34);
    F(F>1e18)=0;
    F=F/rg;
    F(F<1e-2)=0;
    dH = F;
    
    
%    [ZM,ZZ] = sub_zz_zm(fina,finb,HH);
%    ZZ(isnan(ZZ))=100;
%    ZM(isnan(ZM))=100;
    nin=1; % find norm inside the contour, [] - skip this
    A = collocate_UTS_section(II,JJ,U,V,T,S,dH,HH,LON,LAT,nin);
    enorm=sign(A.Norm(:,1)+A.Norm(:,2)); % only sign, as norm comp already defined=un
    un0=A.Unrm; % normal flow
    tn=A.Tnrm;
    sn=A.Snrm;
    DZ=abs(A.dH);
    dx=A.segmL;
    [DX,dmm]=meshgrid(dx,[1:nlr]);

% Normal U:
    clear un
    for kk=1:nlr
      un(kk,:)=un0(kk,:).*enorm';
    end
    
    rhow=sw_dens0(sn,tn);
    
% Adjust bias in overall 
% transport 
% Deep layers have unrealistically high
% velocities (up to 2.5 m/s)
% probably due to time-averaging in HYCOM
% output fields
% This results in unrealistically high
% vol transport that is unbalanced when integrateed 
% along the closed contour
% Adjust deep flow by reducing the high-speed fluxes
    if uadj>0,
      vfmin=-0.1e6; % adjust if total flux < vfmin
      vfmax=1e30; 
      Vadj=sub_crrct_deepU(un,tn,Tref1,Tref2,DX,DZ,Cp,rhow,vfmin,vfmax);
      un=Vadj;
    end
    
    vf0=nansum(un.*DX.*DZ);
    Vflx=nansum(vf0);
    hf1=nansum(un.*Cp.*rhow.*(tn-Tref1).*DX.*DZ);
    Hflx1=nansum(hf1);
    hf2=nansum(un.*Cp.*rhow.*(tn-Tref2).*DX.*DZ);
    Hflx2=nansum(hf2);

    HFLX.VolFlx(cc,:)=vf0;
    HFLX.HFlx1(cc,:)=hf1;
    HFLX.HFlx2(cc,:)=hf2;
    HFLX.TM(cc)=dnmb;

% Box averaged T/S    
    nsc=length(SCT);
    for isc=1:nsc
      IN=SCT(isc).Inpoints;
      dhm=dH(:,IN);
      acl=Acell(IN);
      vcl=[];
      for kk=1:nlr
	vcl(kk,:)=dhm(kk,:).*acl'; % grid cell volume
      end
      vtot=nansum(nansum(vcl));
      tmm=T(:,IN);
      tmn=nansum(nansum(tmm.*vcl))./vtot;
      smm=S(:,IN);
      smn=nansum(nansum(smm.*vcl))./vtot;
      SCT(isc).Tbox(cc)=tmn;
      SCT(isc).Sbox(cc)=smn;
      SCT(isc).TM(cc)=dnmb;
      fprintf('%s Tbox=%6.2f Sbox=%6.2f\n',SCT(isc).Name,tmn,smn);
    end

%    Vmn=Vmn+Vflx;
%    Hmn1=Hmn1+nansum(hf1);
%    Hmn2=Hmn2+nansum(hf2);
    dmm=HFLX.VolFlx;
    Vmn=sum(nanmean(dmm,1));
    dmm=HFLX.HFlx1;
    Hmn1=sum(nanmean(dmm,1));
    dmm=HFLX.HFlx2;
    Hmn2=sum(nanmean(dmm,1));
    fprintf('VFlux=%6.1f Sv, HFlux1=%10.8d, HFlux2=%10.8d\n',...
	    Vflx*1e-6, Hflx1, Hflx2);
    fprintf('Mean %i rec: VFlux=%8.6fSv, HFlx1=%10.8d, HFlx2=%10.8d\n',...
	    cc,Vmn*1e-6, Hmn1, Hmn2);
    
    fprintf('1 day processed %6.3f min\n\n',toc/60);
    
%    keyboard
    
  
  if mod(cc,50)==0;
    fprintf('Saving %s\n',fmat);
    save(fmat,'HFLX','SCT');
  end
  
  end  % day

  fprintf('Saving %s\n',fmat);
  save(fmat,'HFLX','SCT');
  
end   % year





    
    



