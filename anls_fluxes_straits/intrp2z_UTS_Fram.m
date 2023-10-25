% Interpolate to Z levels 
% Fram xsection U,T,S
% Extracted 2D arrays in extract_daily_UTS_Fram.m
%

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

expt = 112;
TV   = 11;  % topo version
sctnm='Fram';

s_mat=1;
YR1=2005;
YR2=2005;

rg  = 9806;
hgg = 1e20; % 

pthmat  = '/Net/mars/ddmitry/hycom/ARCc0.08/data_mat/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
btx='anls_UTS_Fram.m';

btx='inrp2z_UTS_Fram.m';

% Interpolate onto fixed z-grid
% Specify fixed ZZ (interf) and ZM levels:
ZZf = [(0:-1:-50)';(-52:-2:-200)';...
       (-205:-5:-1000)';(-1025:-25:-3050)'];
kzz = length(ZZf);

dZf=diff(ZZf);
ZMf = [];
for ik=1:kzz-1
  ZMf(ik,1)=ZZf(ik)+0.5*(ZZf(ik+1)-ZZf(ik));
end


for yr=YR1:YR2
  fmat_in = sprintf('%sarc08-%3.3i_UTS_daily_%s_%4.4i.mat',...
		 pthmat,expt,sctnm,yr);

  fmat_out = sprintf('%sarc08-%3.3i_UTSonZ_daily_%s_%4.4i.mat',...
		 pthmat,expt,sctnm,yr);
  fprintf('Loading %s\n',fmat_in);
  load(fmat_in);

  SCTZ=struct;

  TM=SCT.TM;
  nrc=length(TM);
  IS=SCT.I;
  JS=SCT.J;
  %X=SCT.X(1:end-1);
  %Y=SCT.Y(1:end-1);
  X=SCT.SGM_Lon;
  Y=SCT.SGM_Lat;
  dL=SCT.SGM_DX;

  f_sct=0;
  if f_sct==1
    ftopo = sprintf('%sdepth_ARCc0.08_%2.2i.nc',pthtopo,TV); % 
    fprintf('Getting topo %s\n',ftopo);
    HH  = nc_varget(ftopo,'Bathymetry');
    LON = nc_varget(ftopo,'Longitude');
    LAT = nc_varget(ftopo,'Latitude');
    [mm,nn]=size(LON);
    [DX,DY]=sub_dx_dy(LON,LAT);

    figure(1); clf;
    contour(HH,[0 0],'k');
    hold on;
    contour(HH,[-5000:1000:-10],'Color',[0.7 0.7 0.7]);
    contour(LAT,[79 79],'b');
    axis('equal');
    set(gca,'xlim',[800 1200],...
	    'ylim',[800 1100]);
    plot(IS,JS,'r.');
    title('Section, 79N');

  end

  SCTZ.ZZ=ZZf;
  SCTZ.ZM=ZMf;
  SCTZ.I=IS;
  SCTZ.J=JS;
  SCTZ.LON=X;
  SCTZ.LAT=Y;

  nrc=length(TM);
  
  for it=1:nrc
    tic;
    dnmb=TM(it);
    dv=datevec(dnmb);
    Hb=SCT.SGM_Hb;
    T=squeeze(SCT.Temp(it,:,:));
    S=squeeze(SCT.Saln(it,:,:));
    V=squeeze(SCT.Unrm(it,:,:));
    dH=squeeze(SCT.dH_thkn(it,:,:));
    [ma,na]=size(dH);
    ZZ=zeros(1,210);
    a=-cumsum(dH,1);
    ZZ=[ZZ;a];
    ll=size(ZZ,1);

    fprintf('%s\n',datestr(dnmb));
    
    ZM=[];
    for kk=1:ll-1
      z1=squeeze(ZZ(kk,:,:));
      z2=squeeze(ZZ(kk+1,:,:));
      ZM(kk,:)=0.5*(z2+z1);
    end

% Check: area-mean interpolated values and original 
% should be the same
    if ~exist('DL','var');
      [mm,nn]=size(V);
      [DL,dmm]=meshgrid(dL,[1:mm]);
    end
    Acell=DL.*dH;
    Ih=find(dH>0);
    Asm=sum(Acell(Ih));
    Vm=1/Asm*sum(V(Ih).*Acell(Ih));
    Tr=sum(V(Ih).*Acell(Ih));
    Sm=1/Asm*sum(S(Ih).*Acell(Ih));
    Tm=1/Asm*sum(T(Ih).*Acell(Ih));
    
    
%    keyboard
  %
  % Interpolate onto z-level 2D 
  % Temperature
    Tz = sub_interp2z_2D(SCT,T,ZZf,ZMf,ZZ,ZM);  
    Tz = sub_add_bottom(Tz,Hb,ZZf);

    Sz = sub_interp2z_2D(SCT,S,ZZf,ZMf,ZZ,ZM);  
    Sz = sub_add_bottom(Sz,Hb,ZZf);

    Vz = sub_interp2z_2D(SCT,V,ZZf,ZMf,ZZ,ZM);  
    Vz = sub_add_bottom(Vz,Hb,ZZf);

    SCTZ.TM(it)=dnmb;
    SCTZ.Temp(it,:,:)=Tz;
    SCTZ.Saln(it,:,:)=Sz;
    SCTZ.NVel(it,:,:)=Vz;
%
% Check if means are preserved
    Azcell=abs(dZf*dL);
    Izh=find(~isnan(Tz));
    Azsm=sum(Azcell(Izh));
    Vzm=1/Azsm*sum(Vz(Izh).*Azcell(Izh));
    Trz=sum(Vz(Izh).*Azcell(Izh));
    Szm=1/Azsm*sum(Sz(Izh).*Azcell(Izh));
    Tzm=1/Azsm*sum(Tz(Izh).*Azcell(Izh));

    errV=abs((Trz-Tr)/Tr);
    errS=abs((Szm-Sm)/Sm);
    errT=abs((Tzm-Tm)/Tm);
    
    fprintf('Check interpolation: Original vs Interpolated, error:\n');
    fprintf('Mean transp, Sv: %8.2f, %8.2f %6.4d\n', Tr*1e-6,Trz*1e-6,errV);
    fprintf('Mean salin,    : %8.2f, %8.2f %6.4d\n', Sm,Szm,errS);
    fprintf('Mean temp,  dgr: %8.2f, %8.2f %6.4d\n', Tm,Tzm,errT);
    
    if (errV>0.1),
%      error('Transport not preserved ...');
      fprintf(' !!!! Transport not preserved !!! \n\n');
    end
    
    f_chck=0;
    if f_chck==1
      nfg=2;
      [ma,na]=size(dH);

      stl=sprintf('arc08-%3.3i %s, T, %4.4i/%2.2i/%2.2i',...
		  expt,sctnm,dv(1:3));
      cntr=[-2:0.5:10];
      cc0=0;
      mbtm=0;
      sub_plot_Txsct(nfg,Hb,T,ZZ,X,dnmb,stl,cntr,cc0,mbtm);
      bottom_text(btx,'pwd',1);


      stl=sprintf('arc08-%3.3i %s, T on Z, %4.4i/%2.2i/%2.2i',...
		  expt,sctnm,dv(1:3));
      nfg=3;
      mbtm=0;  % add bottom patch
      sub_plot_Txsct(nfg,Hb,Tz,ZMf,X,dnmb,stl,cntr,cc0,mbtm);
      bottom_text(btx,'pwd',1);


      stl=sprintf('arc08-%3.3i %s, S, %4.4i/%2.2i/%2.2i',...
		  expt,sctnm,dv(1:3));
      cntr=[31:0.5:35];
      cc0=34.8;
      nfg=4;
      mbtm=0;
      sub_plot_Sxsct(nfg,Hb,S,ZZ,X,dnmb,stl,cntr,cc0,mbtm);
      bottom_text(btx,'pwd',1);

      stl=sprintf('arc08-%3.3i %s, S on Z, %4.4i/%2.2i/%2.2i',...
		  expt,sctnm,dv(1:3));
      nfg=5;
      mbtm=0;
      sub_plot_Sxsct(nfg,Hb,Sz,ZMf,X,dnmb,stl,cntr,cc0,mbtm);
      bottom_text(btx,'pwd',1);


    % V component
    % Project on a 79 Lat to avoid "steps" discontinuities
    % due to zig-zaging segments
      dI=diff(SCT.I);
      I0=find(dI>0);
      V0=V(:,I0);
      V0z=Vz(:,I0);
      Hb0=Hb(I0);
      ZZ0=ZZ(:,I0);
      X0=X(I0);

      stl=sprintf('arc08-%3.3i %s, V, %4.4i/%2.2i/%2.2i',...
		  expt,sctnm,dv(1:3));
      cntr=[-0.5:0.1:0.5];
      cc0=0;
      nfg=6;
      mbtbm=0;
      sub_plot_Vxsct(nfg,Hb0,V0,ZZ0,X0,dnmb,stl,cntr,cc0,mbtbm);
      bottom_text(btx,'pwd',1);

	stl=sprintf('arc08-%3.3i %s, V, %4.4i/%2.2i/%2.2i',...
		  expt,sctnm,dv(1:3));
      nfg=7;
      mbtbm=0;
      sub_plot_Vxsct(nfg,Hb0,V0z,ZMf,X0,dnmb,stl,cntr,cc0,mbtbm);
      bottom_text(btx,'pwd',1);
    end

    fprintf('Interpolated 1 rec, %7.3f sec\n',toc);
  %keyboard  
  end

  if s_mat==1
    fprintf('Saving %s\n',fmat_out);
    save(fmat_out,'SCTZ');
  end

end  % YR






