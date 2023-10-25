% river passive tracers from ARCc0.04
% distributed along the Greenland Coast
% Not depth-integrated 
% Analyze amount of tracers along the coast
% In two models - Trying to understand
% why in the coarser model there is more 
% tracer advected away off the coast?
% Hypothesis is - in higher resolution model
% The trapping is much stronger and more
% tracer is kept along the coast
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

plr   = 1;  % layer to plot
s_mat = 0;
f_extr_sgm = 0;
s_reg = 0; %=0 - archm files are whole domain
         %=1 - archm files subsampled into nAtl. Region,
	 %      check fortran/hycom/extract_subdomain.F90
nTr=1;   % tracer to plot
TrMn=1e-10; % Threshold value to plot
rg = 9806;
fprintf('Tracer #: %i, Threhold value: %8.5d\n',nTr,TrMn);


YRPLT=[];
cc=0;
for iyr=2005:2005
  for idd=210:210
%for iyr=2008:2008
%  for idd=1:7:365
%    if idd==1, idd=2; end;
    cc=cc+1;
    YRPLT(cc,1)=iyr;
    YRPLT(cc,2)=idd;
  end
end

np=size(YRPLT,1);

% Experiments:
% 102 - test simulation, HYCOM GLBb0.08 nest: 1 file kept constant (1993,1,1)

regn = 'ARCc0.04';
expt = 011;  
%regn = 'ARCc0.08';
%expt = 110;  

pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthmat = '/Net/ocean/ddmitry/HYCOM/ARCc/data_mat/';
%pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/topo_grid/';


inc1=1;
inc2=3200;
jnc1=1;
jnc2=5040;
djnc=5040;
dinc=3200;

%inc1=1;
%inc2=1600;
%jnc1=1;
%jnc2=2520;
%djnc=2520;
%dinc=1600;


ftopo = sprintf('%s/depth_%s_17DD.nc',pthtopo,regn); % 
%ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');

[mm,nn]=size(LON);

%hmsk=HH;
%hmsk(HH<0)=nan;

[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2

xlim1 = 0;
xlim2 = nn;
ylim1 = 0;
ylim2 = mm;

fmt04 = sprintf('%sarc04_green_segm.mat',pthmat);
if f_extr_sgm==1
% Read in segments from ARCc0.08
% see trcr_GreenCoast.m
  fmt08 = sprintf('%sarc08_green_segm.mat',pthmat);
  fprintf('Loading %s\n',fmt08);
  load(fmt08);

  X8=SGM.X;
  Y8=SGM.Y;
  Xo8=SGM.Xout;
  Yo8=SGM.Yout;
  nx=length(X8);
  fprintf('Finding coordinates of Green. contours ...\n');
  clear IIs JJs Iout Jout
  for ik=1:nx
    if mod(ik,100)==0;
      fprintf('#1  ... Done %6.3f %%\n',(ik/nx)*100);
    end

    x=X8(ik);
    y=Y8(ik);
    dmm=distance_spheric_coord(y,x,LAT,LON);
    dmm(dmm>10000)=nan;
    I=find(dmm==min(min(dmm)),1);
    [jj,ii] = ind2sub(size(HH),I);
    IIs(ik) = ii;
    JJs(ik) = jj;
  end;
  nx = length(Xo8);
  for ik=1:nx
    if mod(ik,100)==0;
      fprintf('#2  ... Done %6.3f %%\n',(ik/nx)*100);
    end
    x=Xo8(ik);
    y=Yo8(ik);
    dmm=distance_spheric_coord(y,x,LAT,LON);
    dmm(dmm>10000)=nan;
    I=find(dmm==min(min(dmm)),1);
    [jj,ii] = ind2sub(size(HH),I);
    Iout(ik) = ii;
    Jout(ik) = jj;
  end
  clear SGM;
  SGM.I=IIs;
  SGM.J=JJs;
  SGM.Iout=Iout;
  SGM.Jout=Jout;
  fmt04 = sprintf('%sarc04_green_segm.mat',pthmat);
  fprintf('Saving %s\n',fmt04);
  save(fmt04,'SGM');
else
  fprintf('Loading %s\n',fmt04);
  load(fmt04);
end  
  
IIs = SGM.I;
JJs = SGM.J;
IIc = SGM.Iout;
JJc = SGM.Jout;
INDs= sub2ind(size(HH),JJs,IIs);
INDc= sub2ind(size(HH),JJc,IIc);


% Check segments:
f_chck=0;
if f_chck==1
  figure(1); clf;
  contour(HH,[0 0],'k');
  hold on;
  contour(HH,[-500 -500],'g');
  plot(IIs,JJs,'r.')
  plot(IIc,JJc,'.');
end
%keyboard


mdx = mean(DX(INDs));
mdy = mean(DY(INDs));
dxy = sqrt(mdx^2+mdy^2);
dstm = 300e3; % distance offshore 
dst = round(dstm/dxy); % aprx # of grid points to get dst length

% Read fields:
cnc=0;
ip1=1;
for ip=ip1:np
  yr=YRPLT(ip,1);
  iday=YRPLT(ip,2);
  pthbin = sprintf('/nexsan/hycom/ARCc0.04_011/data/%i/',yr);  
%  pthbin = sprintf('/nexsan/archive/%s_%3.3i/data/%i/',regn,expt,yr)

  fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
  finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);

%  fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
%  finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);
  
  if ~exist(fina,'file');
    fprintf('Not found: %s\n\n',fina);
    continue;
  end
  
  cnc=cnc+1;
  dnmb=datenum(yr,1,1)+iday-1;
  DV=datevec(dnmb);

%fprintf('Plotting: Tracer #%i V.Layer=%i\n\n',nTr,plr);
  fprintf('%4.4i_%2.2i_%2.2i: %s\n',DV(1:3),fina);
% in archm u_vel=utot; total velocity
%  keyboard
  f_getu=0;
  if f_getu==1
    ilr=1;
    [Fu,n,m,l] = read_hycom(fina,finb,'u-vel.');
    [Fv,n,m,l] = read_hycom(fina,finb,'v-vel.');
    U=squeeze(Fu(ilr,:,:));
    U(U>1e6)=nan;
    V=squeeze(Fv(ilr,:,:));
    V(V>1e6)=nan;
    S=sqrt(U.^2+V.^2);
  
  end

 % Layer thickness:
  [F,n,m,l] = read_hycom(fina,finb,'thknss','r_layer',plr);
%  F=squeeze(F(plr,:,:));
  F=F./rg;
  F(F>1e10)=0;
  F(F<1e-2)=0;
  dH=squeeze(F); 

  tic;
  [F,n,m,l] = read_hycom(fina,finb,'tracer','r_tracer',nTr,'r_layer',plr);
  toc;
  F(F>1e6)=nan;

  if plr==1
    Tr=squeeze(F(plr,:,:));   % tr. conc, kg/m3
    mTr=squeeze(F(plr,:,:)).*dH.*DX.*DY; % kg - tracer mass in 1 grid cell
  else 
    smm = HH*0;
    shh = HH*0;
    for kk=1:plr
      tr = squeeze(F(kk,:,:));
      dh = squeeze(dH(kk,:,:));
      smm = smm+tr.*dh;
      shh = shh+dh;
    end
    mTr = smm.*DX.*DY; % mass, kg, depth-cell intergrated
    Tr = tr./shh;   % tr. conc, kg/m3, depth-averaged
  end
    
% Threshold value:
  Tr(Tr<=TrMn)=nan;
  I=find(~isnan(Tr));
  
  nlr = 1; % # of layers to average
  clear i
%  rtr = i; % where to turn to get off land
  dd = 5;  % points to skip along coast segments
  TRC = sub_get_coastTr(HH,dH,Tr,mTr,LON,LAT,SGM,plr,dst,dd,dnmb);
  
end;  % day loop

%DL = TRC.Dist_m;
%sTr = TRC.Tr;
if s_mat==1
  fmt = sprintf('%s%s_trcr_NLr%2.2i_GreenCoast_%i.mat',pthmat,regn,plr,yr);
  fprintf('Saving %s\n',fmt);
  save(fmt,'TRC');
end






