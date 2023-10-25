% Interpolate HYCOM 0.04 fields onto
% AWI grid
% for calculating vorticity etc
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

expt = '011';
TV = '17DD';

YR1 = 2005;
YR2 = 2007;

rg=9806;  % convert pressure to depth, m
huge=1e10;
omg = 7.2921159e-5; 


sfig = 0;
zz0  = -100;  %  depth for interpolation

pthmat  = '/Net/ocean/ddmitry/HYCOM/ARCc/data_eddy_AWI/';
pthout  = '/Net/tholia/ddmitry/hycom/ARCc0.04/data_awi_intrp/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/topo_grid/';
pthfig  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/%s/fig_2D/',...
		  expt);

% Get AWI grid:
fmat = sprintf('%sgrid_eddy_tracking.mat',pthmat);
fprintf('Loading %s\n',fmat);
AWI = load(fmat);
alat=AWI.lat;
alon=AWI.lon;
lt1=min(alat);
lt2=max(alat);
ln1=min(alon);
ln2=max(alon);


% Barycentric interpolation weights and HYCOM indices
fout    = sprintf('%sindx2awi004.mat',pthmat);
fprintf('Loading %s\n',fout);
load(fout);


ftopo = sprintf('%sdepth_ARCc0.04_%s.nc',pthtopo,TV); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

% Find AWI domain
dd=distance_spheric_coord(LAT,LON,lt1,ln1);
[j11,i11]=find(dd==min(min(dd)));
dd=distance_spheric_coord(LAT,LON,lt2,ln1);
[j21,i21]=find(dd==min(min(dd)));
dd=distance_spheric_coord(LAT,LON,lt1,ln2);
[j12,i12]=find(dd==min(min(dd)));
dd=distance_spheric_coord(LAT,LON,lt2,ln2);
[j22,i22]=find(dd==min(min(dd)));
JH1=min([j11,j12,j21,j22])-2;
JH2=max([j11,j12,j21,j22])+2;
IH1=min([i11,i12,i21,i22])-2;
IH2=max([i11,i12,i21,i22])+2;


f_map=0;
if f_map==1
  figure(1); clf;
  contour(HH,[0 0],'k'); hold
  contour(HH,[zz0 zz0],'r');
  contour(LON,[-170:10:170],'Color',[0.5 0.5 0.5]);
  contour(LAT,[40:10:79],'Color',[0.5 0.5 0.5]);
  axis('equal');
  lt1=min(alat);
  lt2=max(alat);
  ln1=min(alon);
  ln2=max(alon);

  contour(LON,[ln1 ln1],'g');
  contour(LON,[ln2 ln2],'g');
  contour(LAT,[lt1 lt1],'g');
  contour(LAT,[lt2 lt2],'g');
  title('AWI Fram domain');
  
      
end

f_chck=0;
if f_chck==1
  i0=100;
  j0=80;
  x0=alon(i0);
  y0=alat(j0);
  im1=INDX(1).hycom_I(j0,i0);
  jm1=INDX(1).hycom_J(j0,i0);
  im2=INDX(2).hycom_I(j0,i0);
  jm2=INDX(2).hycom_J(j0,i0);
  im3=INDX(3).hycom_I(j0,i0);
  jm3=INDX(3).hycom_J(j0,i0);
  
  lmb1=INDX(1).brctr_lmbd(j0,i0);
  lmb2=INDX(2).brctr_lmbd(j0,i0);
  lmb3=INDX(3).brctr_lmbd(j0,i0);
  
  figure(10); clf;
  plot(x0,y0,'r*');
  hold on;
  x1=LON(jm1,im1);
  x2=LON(jm2,im2);
  x3=LON(jm3,im3);
  y1=LAT(jm1,im1);
  y2=LAT(jm2,im2);
  y3=LAT(jm3,im3);
  plot([x1 x2],[y1 y2],'b-');
  plot([x1 x3],[y1 y3],'b-');
  plot([x3 x2],[y3 y2],'b-');
  text(x1,y1,sprintf('%4.3f',lmb1));
  text(x2,y2,sprintf('%4.3f',lmb2));
  text(x3,y3,sprintf('%4.3f',lmb3));
  
end

[DX,DY]=sub_dx_dy(LON,LAT);
Fcor = 2*omg*sind(LAT);
Lmsk = HH*0;
Lmsk(HH<zz0) = 1;

for yr=YR1:YR2
  mold=0;
  cc=0;
%  iday1 = datenum(yr,12,1)-datenum(yr,1,1)+1;
  iday1 =1;
  iday2 = datenum(yr+1,1,1)-datenum(yr,1,1);
%  iday2 = 31;
  for iday=iday1:iday2
    tic;
    dnmb=datenum(yr,1,1)+iday-1;
    DV =datevec(dnmb);
    yr = DV(1);
    mo = DV(2);
    if mold==0, mold=mo; end;
    
    fprintf('\n---   %4.4i/%2.2i/%2.2i ---\n',DV(1:3));

    pthbin = sprintf('/nexsan/hycom/ARCc0.04_%s/data/%4.4i/',expt,yr);

    fina = sprintf('%s%s_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
    finb = sprintf('%s%s_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);

    if ~exist(fina,'file') | ~exist(finb,'file');
      fprintf('Not found *a or *b: %s\n',fina);
      fprintf('                     %s\n',finb);
     continue;
    end

    % find layer  with zz0
    if ~exist('iz0','var')
      [ZM,ZZ] = sub_zz_zm(fina,finb,HH);
      ll=size(ZM,1);
      ZZ(isnan(ZZ))=100;
      ZM(isnan(ZM))=100;
      A=ZM(:,JH1:JH2,IH1:IH2);
      hs = HH(JH1:JH2,IH1:IH2);
      hs(hs>zz0)   = nan;
      A(isnan(hs)) = nan;
      A(A>0)       = nan;
      zlr = [];
      for k=1:ll
        dmm=squeeze(A(k,:,:));
        zlr(k,1)=nanmean(nanmean(dmm));
      end
      dmm=abs(zlr-zz0);
      iz0=find(dmm==min(dmm));
      fprintf('Closest layer %i = %4.1fm to %4.1fm\n',iz0,zlr(iz0),zz0);
    end      

% Get pang to rotate U,V components to a new grid:
    if ~exist('pang','var')
      fgrda = sprintf('%sregional.grid.a',pthtopo);
      fgrdb = sprintf('%sregional.grid.b',pthtopo);
      [pang,a1,a2] = read_pang(fgrda,fgrdb);
    end
    
    
    cc = cc+1;
    fprintf('Getting data expt %s: %4.4i_%2.2i_%2.2i: %s\n',expt,DV(1:3),fina);

    [F,n,m,nlr] = read_hycom(fina,finb,'u-vel.','r_layer',iz0);
    F(F>huge)=nan;
    U=squeeze(F);

    [F,n,m,nlr] = read_hycom(fina,finb,'v-vel.','r_layer',iz0);
    F(F>huge)=nan;
    V=squeeze(F);

    [F,n,m,nlr] = read_hycom(fina,finb,'temp','r_layer',iz0);
    F(F>huge)=nan;
    T=squeeze(F);
 
    [F,n,m,nlr] = read_hycom(fina,finb,'salin','r_layer',iz0);
    F(F>huge)=nan;
    S=squeeze(F);
    
% Rotate from X, Y hycom grid to 
% North/East components
    Ur = cos(-pang).*U + sin(-pang).*V;
    Vr = cos(-pang).*V - sin(-pang).*U;
    
    f_chu=0;
    if f_chu==1
      i=1550;
      j=2200;
      uu=U(j,i);
      vv=V(j,i);
      ur=Ur(j,i);
      vr=Vr(j,i);
      
      figure(11); clf;
      compass(uu,vv,'r');
      hold on
      compass(ur,vr,'b'); % rotated component
     end
     clear U V   
      
      
%    keyboard
% Interpolation into AWI grid:
     im1 = INDX(1).hycom_I; % i HYCOM index of vertex 1
     jm1 = INDX(1).hycom_J; % j HYCOM index of vertex 1
     im2 = INDX(2).hycom_I;
     jm2 = INDX(2).hycom_J;
     im3 = INDX(3).hycom_I;
     jm3 = INDX(3).hycom_J;
     [ma,na] = size(im1);
     im1 = reshape(im1,ma*na,1);
     jm1 = reshape(jm1,ma*na,1);
     im2 = reshape(im2,ma*na,1);
     jm2 = reshape(jm2,ma*na,1);
     im3 = reshape(im3,ma*na,1);
     jm3 = reshape(jm3,ma*na,1);
     I1  = sub2ind([mm,nn],jm1,im1);
     I2  = sub2ind([mm,nn],jm1,im1);
     I3  = sub2ind([mm,nn],jm1,im1);
     
     lmb1 = INDX(1).brctr_lmbd;
     lmb2 = INDX(2).brctr_lmbd;
     lmb3 = INDX(3).brctr_lmbd;
     lmb1 = reshape(lmb1,ma*na,1);
     lmb2 = reshape(lmb2,ma*na,1);
     lmb3 = reshape(lmb3,ma*na,1);

     Ui = lmb1.*Ur(I1)+lmb2.*Ur(I2)+lmb3.*Ur(I3);
     Vi = lmb1.*Vr(I1)+lmb2.*Vr(I2)+lmb3.*Vr(I3);
     Ti = lmb1.*T(I1)+lmb2.*T(I2)+lmb3.*T(I3);
     Si = lmb1.*S(I1)+lmb2.*S(I2)+lmb3.*S(I3);

     Ui = reshape(Ui,[ma,na]);
     Vi = reshape(Vi,[ma,na]);
     Ti = reshape(Ti,[ma,na]);
     Si = reshape(Si,[ma,na]);

     HYCOM.Time(cc)=dnmb;
     HYCOM.U(cc,:,:)=Ui;
     HYCOM.V(cc,:,:)=Vi;
     HYCOM.T(cc,:,:)=Ti;
     HYCOM.S(cc,:,:)=Si;

     f_chip = 0;
     if f_chip==1
%       A = Ur(JH1:JH2,IH1:IH2);
       A = Vr(JH1:JH2,IH1:IH2);
       [mk,nk] = size(A);
       figure(15); clf;
       subplot(2,1,1);
       pcolor(A); shading flat; % HYCOM field
       colorbar
       axis('equal');
%       title('HYCOM U, rotated to N/E AWI grid');
       title('HYCOM V, rotated to N/E AWI grid');
       caxis([-0.4 0.4]);
       set(gca,'xlim',[0 nk],...
	       'ylim',[0 mk]);
       
       subplot(2,1,2);
%       pcolor(Ui); shading flat;
       pcolor(Vi); shading flat;
       colorbar
       axis('equal');
%       title('U Interp  AWIgrid');
       title('V Interp  AWIgrid');
       caxis([-0.4 0.4]);
       set(gca,'xlim',[0 na],...
	       'ylim',[0 ma]);
      
       btx = 'interp2awi004.m';
       bottom_text(btx,'pwd',1);
       
       
       A = S(JH1:JH2,IH1:IH2);
       [mk,nk] = size(A);
       figure(16); clf;
       subplot(2,1,1);
       pcolor(A); shading flat; % HYCOM field
       colorbar
       axis('equal');
       title('HYCOM S, rotated to N/E AWI grid');
%       caxis([-2 8]);
       caxis([31 35]);
       set(gca,'xlim',[0 nk],...
	       'ylim',[0 mk]);
       
       subplot(2,1,2);
%       pcolor(Ti); shading flat;
       pcolor(Si); shading flat;
       colorbar
       axis('equal');
       title('S Interp  AWIgrid');
%       caxis([-2 8]);
       caxis([31 35]);
       set(gca,'xlim',[0 na],...
	       'ylim',[0 ma]);
      
       btx = 'interp2awi004.m';
       bottom_text(btx,'pwd',1);
       
     end
       
     fprintf('1 record processed %6.1f min\n\n',toc/60);

% Check next date for saving output:
    dnxt = dnmb+1;
    dvn  = datevec(dnxt);
    if mold~=dvn(2)
      fintrp = sprintf('%sarc004_%s_archm2awi.%4.4i_%2.2i.mat',...
		       pthout,expt,yr,mo);
      fprintf('Saving %s\n\n',fintrp);
      save(fintrp,'HYCOM');
      clear HYCOM
      HYCOM = struct;
      cc=0;
      mold=0;
    end


  end
end

if cc>0 & yr==YR2
  fintrp = sprintf('%sarc004_%s_archm2awi.%4.4i_%2.2i.mat',...
		   pthout,expt,yr,mo);
  fprintf('Saving %s\n\n',fintrp);
  save(fintrp,'HYCOM');
  cc=0;
  mold=0;
end

