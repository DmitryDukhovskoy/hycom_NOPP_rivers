% Plot sections
% analyze fluxes:
% Plot seasonal(or monthly) sections of U,T,S
% and  Vol, heat, FW fluxes
%
% averaged over specified years
% across sections on the Gr Shelf
% extracted in flux_Grshelf_sections.m
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

YR1=1993;
YR2=2016;
YRS = [YR1:YR2];
  
Sref=34.9; % Similar to Le Bras and Straneo et al 2018

% =============
% Plot/anls fields:
% =============
f_pltuv=0; % U normal to the section
f_plts=1;  % S
f_pltt=0;  % T
f_fwf=0;   % calc FW flux. mSv
f_vflx=0;   % calc vol flux. Sv
f_uvctr=150; % plot mean U vectors, upper N meters = f_uvctr>0
s_fig=0; 


Tref=0;
rg = 9806;
rhow=1027; 
hgg=1e20;

regn = 'ARCc0.08';
%expt = 110;  % no Gr FWF
expt = 112;  % Gr FWF
pthfig  = '/nexsan/people/ddmitry/Net_ocean/hycom/ARCc0.08/110/fig_GrSct/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat =sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
pthmat2=sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_GrSect/',expt);
btx = 'flux_anls_Grshelf_sections.m';

fprintf('Plotting %s-%3.3i Monthly Greenl Shelf Sections %i-%i\n\n',...
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




fmatu=sprintf('%sarc08_expt%3.3i_greenl_shsect_mnthUzlv_%i.mat',pthmat,expt,2016);
load(fmatu);
ZZ=UZGR(1).ZZlevels;
ZM=UZGR(1).ZM;
nzm=length(ZM);
nzz=length(ZZ);
dZ=abs(diff(ZZ));

f_map=0;
if f_map==1
  fn=5;
  sub_plot_Greenl_contour(HH,LON,LAT,fn,GC);

  for ip=1:nsct
%    IJ=SCT(ip).IJ;
%    plot([IJ(1,1) IJ(2,1)],[IJ(1,2) IJ(2,2)],...
%	 'Linewidth',2.5,'Color',[1. 0.6 0]);
    IIs=UZGR(ip).I;
    JJs=UZGR(ip).J;
    plot(IIs,JJs,'-',...
	 'Linewidth',2.5,'Color',[1. 0.6 0]);
  end
  bottom_text(btx,'pwd',1);
end



cc=0;
for YR=YR1:YR2
  yr=YR;
  fmatu=sprintf('%sarc08_expt%3.3i_greenl_shsect_mnthUzlv_%i.mat',pthmat,expt,YR);
  fmatt=sprintf('%sarc08_expt%3.3i_greenl_shsect_mnthTzlv_%i.mat',pthmat,expt,YR);
  fmats=sprintf('%sarc08_expt%3.3i_greenl_shsect_mnthSzlv_%i.mat',pthmat,expt,YR);
  fmatuv=sprintf('%sarc08_expt%3.3i_greenl_shsect_mnthUVz_%i.mat',pthmat,expt,YR);


  
% Volume flux:
  if ~exist(fmatu,'file')
    fprintf('Year not found: %i, %s\n',YR,fmatu);
    continue;
  end
  
  fprintf('Loading %s\n',fmatu);
  load(fmatu);

  if ~exist('nsct','var') 
    nsct=length('UZGR');
  end      
  
  cc=cc+1;
%
% UZGR is Normal velocity (normal to small segments 
% - only useful for flux calculation)
% UV - U,V components along section
  for ip=1:nsct
    if cc==1;
      SM(ip).usm=UZGR(ip).U;
    else
      SM(ip).usm=SM(ip).usm+UZGR(ip).U;
    end
  end
  
% U,V components:
  fprintf('Loading %s\n',fmatuv);
  load(fmatuv);
  for ip=1:nsct
    if cc==1;
      UV(ip).usm=UVZGR(ip).Ucomp;
      UV(ip).vsm=UVZGR(ip).Vcomp;
    else
      UV(ip).usm=UV(ip).usm+UVZGR(ip).Ucomp;
      UV(ip).vsm=UV(ip).vsm+UVZGR(ip).Vcomp;
    end
  end

  
% S
  fprintf('Loading %s\n',fmats);
  load(fmats);
  for ip=1:nsct
    if cc==1;
      SZ(ip).ssm=SZGR(ip).S;
    else
      SZ(ip).ssm=SZ(ip).ssm+SZGR(ip).S;
    end
  end

% T
  fprintf('Loading %s\n',fmatt);
  load(fmatt);
  for ip=1:nsct
    if cc==1;
      SZ(ip).tsm=TZGR(ip).T;
    else
      SZ(ip).tsm=SZ(ip).tsm+TZGR(ip).T;
    end
  end

% Calculate FW flux
  for ip=1:nsct
    dL=UZGR(ip).Segm_dL;
    I=UZGR(ip).I;
    J=UZGR(ip).J;

    iGr=UZGR(ip).GrCntr_I;
    jGr=UZGR(ip).GrCntr_J;

    dxx=sqrt((I-iGr).^2+(J-jGr).^2);
    Ix=find(dxx==0);
  %  Ix=round(Ix*0.6); - roughly where CF2 station GrShelf in Bras & Straneo 2018

    [DX,DZ]=meshgrid(dL(1:Ix),dZ);
    [DXt,DZt]=meshgrid(dL,dZ);

    for im=1:12
      ss=squeeze(SZGR(ip).S(im,:,1:Ix)); % S 12mo x Vert levels x Segments
      un=squeeze(UZGR(ip).U(im,:,1:Ix)); % U normal
      fwf=((Sref-ss)./Sref).*un.*DX.*DZ;
      FWF(ip).FWF_GrSh(cc,im)=nansum(nansum(fwf)); %only Gr Shelf - contour
      vtr=un.*DX.*DZ;
      VTR(ip).VTR_GrSh(cc,im)=nansum(nansum(vtr)); %only Gr Shelf - contour
% Whole section
      ss=squeeze(SZGR(ip).S(im,:,:)); % S 12mo x Vert levels x Segments
      un=squeeze(UZGR(ip).U(im,:,:)); % U normal
      fwf=((Sref-ss)./Sref).*un.*DXt.*DZt;
      FWF(ip).FWF_sect(cc,im)=nansum(nansum(fwf)); % whole section
      vtr=un.*DXt.*DZt;
      VTR(ip).VTR_sect(cc,im)=nansum(nansum(vtr));

    end    
%    keyboard
  end

% T
  
  
end;   % years    

for ip=1:nsct
  Un=SM(ip).usm/cc;
  SM(ip).usm=Un;
  
  Un=UV(ip).usm/cc;
  UV(ip).usm=Un;
  
  Un=UV(ip).vsm/cc;
  UV(ip).vsm=Un;
  
  Un=SZ(ip).ssm/cc;
  SZ(ip).ssm=Un;
  
  Un=SZ(ip).tsm/cc;
  SZ(ip).tsm=Un;
  
end

% Define winter and summer :
js1=7;
js2=9;
jw1=10;
jw2=3;

% Project U on normal plain to get rid of
% "steps" in +/- flux across zigzagin transect
for ik=1:4
  LL=UZGR(ik).Dist_origin*1e-3; % m
  dx=3.0;
  dL=diff(LL);
  dL=[dL;dL(end)];
  Ip=find(dL>dx);
  LLu=LL(Ip);
  
  dmm=SM(ik).usm;
  Us=squeeze(nanmean(dmm(js1:js2,:,:),1));
  Us=Us(:,Ip);
  
  I=[1:jw1,jw2:12];
  Uw=squeeze(nanmean(dmm(I,:,:),1));
  Uw=Uw(:,Ip);
 
  SM(ik).Usummer=Us;
  SM(ik).Uwinter=Uw;
end

 
 

% For plotting add an extra layer 
%Un(nzz,:)=Un(nzz-1,:);

dw=0.4;
dh=0.38;
POS=[0.07 0.56 dw dh; ...
     0.53 0.56 dw dh; ...
     0.07 0.08 dw dh; ...
     0.53 0.08 dw dh];




% =================
% Plotting UV
% =================
if f_pltuv==1
  stlS=sprintf('arc08-%3.3i, U m/s, %i-%i, mo:%i-%i',expt,YR1,YR2,js1,js2);
  stlW=sprintf('arc08-%3.3i, U m/s, %i-%i, mo:%i-%i',expt,YR1,YR2,jw1,jw2);
  sub_plot_UV_GrShelf(nf,UV,UZGR,js1,js2,jw1,jw2,btx,stlS,stlW,POS);
end
 

% ============================
% Salinity
% ============================
if f_plts==1
% Summer and winter months
%  js1=7;
%  js2=9;
%  jw1=11;
%  jw2=3;

  nf=3;
  stlS=sprintf('008-%3.3i, S, %i-%i, mo:%i-%i',expt,YR1,YR2,js1,js2);
  stlW=sprintf('008-%3.3i, S, %i-%i, mo:%i-%i',expt,YR1,YR2,jw1,jw2);
  sub_plot_S_GrShelf(nf,SZ,SZGR,js1,js2,jw1,jw2,btx,stlS,stlW,POS);
end

% ============================
% Temperature
% ============================
if f_pltt==1
% Summer and winter months
%  js1=7;
%  js2=9;
%  jw1=11;
%  jw2=3;

  nf=7;
  stlS=sprintf('008-%3.3i, T, %i-%i, mo:%i-%i',expt,YR1,YR2,js1,js2);
  stlW=sprintf('008-%3.3i, T, %i-%i, mo:%i-%i',expt,YR1,YR2,jw1,jw2);
  sub_plot_T_GrShelf(nf,SZ,TZGR,js1,js2,jw1,jw2,btx,stlS,stlW,POS);
end




% ===============
% Plot FW flux
% ===============
if f_fwf==1

  nf=5;
  stl=sprintf('arc08-%3.3i FWF mSv, S=%3.1f',expt,Sref);
  sub_FWF_GrShelf(FWF,UZGR,nf,stl);
  bottom_text(btx,'pwd',1);
  
  if s_fig==1
    fgnm=sprintf('%sarc08_%3.3i_FWF_flux_GrShSect_S34p9',...
	    pthfig,expt);
%    fgnm=sprintf('%sarc08_%3.3i_FWF_flux_GrShSect_S34p8',...
%	    pthfig,expt);
    fprintf('Saving %s\n',fgnm);
    print('-depsc2',fgnm);
  end
  
end
% ===============
% Plot Vol flux
% ===============
if f_vflx==1

  nf=7;
  stl=sprintf('arc08-%3.3i VFlx Sv',expt);
  sub_VolFlx_GrShelf(VTR,UZGR,nf,stl);
  bottom_text(btx,'pwd',1);
  
  if s_fig==1
    fgnm=sprintf('%sarc08_%3.3i_Vol_flux_GrShSections',...
	    pthfig,expt);
    fprintf('Saving %s\n',fgnm);
    print('-depsc2',fgnm);
  end
  
end
