% 0.04 HYCOM-CICE
%
% Plot, Extract T,S,U 
% Daily values
% across any section along X or Y model axis
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

dday=3; 
YR1=2005;
YR2=2005;
id1=241;
id2=241;
YRS = [YR1:YR2];

s_fig=1; 
sctnm='mvp1'; % Amundsen, MVP line 1
%sctnm='mvp2'; % Amundsen, MVP line 2

rg = 9806;
rhow=1027; 
hgg=1e20;

regn = 'ARCc0.04';
%expt = 011;
expt = 012; % Greenland runoff
pthfig  = '/Net/mars/ddmitry/hycom/ARCc0.04/012/fig_xsctTS/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/topo_grid/';
%pthmat =sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
%pthmat = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/%3.3i/data_GrSect/',expt);
btx='TSU_xsct004.m';

fprintf('%s-%3.3i Daily Greenl Shelf SE Section U,T,S Zlevels, %i-%i\n\n',...
	regn,expt,YR1,YR2);

ftopo = sprintf('%s/depth_%s_17DD.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);

% Amundsen MVP section line 1:
switch (sctnm),
 case('mvp1'); % these are exact coord, but will extent the section a bit later
  ln2=-(57+44/60);
  lt2=57+42/60;
  ln1=-(59+32/60);
  lt1=57+56/60;
 case('mvp2');
% Amundsen MVP section line 1:
  ln1=-(57+44/60);
  lt1=57+42/60;
  ln2=-(59+32/60);
  lt2=57+56/60;
end

SCT.Name=sctnm;
IJ=sub_XY2indx([ln1,lt1],LON,LAT);
SCT.IJ(1,:)=IJ;
IJ=sub_XY2indx([ln2,lt2],LON,LAT);
SCT.IJ(2,:)=IJ;

%  For MVP: extend sections and make along X axis
SCT.IJ(1,1)=770;
SCT.IJ(2,1)=900;
SCT.IJ(2,2)=SCT.IJ(1,2);

fprintf('Section: %s\n',SCT.Name);
IJs=SCT.IJ;
[IIs,JJs]=sub_xsct_indx(IJs(1,1),IJs(1,2),IJs(2,1),IJs(2,2));
SCT.I=IIs;
SCT.J=JJs;
%keyboard
nsct=length(SCT);

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
  figure(10); clf;
  contour(HH,[-5000:500:0],'k');
  hold on;

  for ip=1:nsct
%    IJ=SCT(ip).IJ;
%    plot([IJ(1,1) IJ(2,1)],[IJ(1,2) IJ(2,2)],...
%	 'Linewidth',2.5,'Color',[1. 0.6 0]);
    IIs=SCT(ip).I;
    JJs=SCT(ip).J;
    plot(IIs,JJs,'-',...
	 'Linewidth',2.5,'Color',[0.8 0.3 0]);
  end
  axis('equal');
  set(gca,'xlim',[350 950],...
          'ylim',[600 1200]);
  bottom_text(btx,'pwd',1);
keyboard
end

% Interpolate onto fixed z-grid
% Specify fixed ZZ (interf) and ZM levels:
ZZf = [(0:-1:-200)';...
       (-202:-2:-1000)';(-2025:-25:-3000)';(-3050:-50:-5000)'];
kzz = length(ZZf);

dZf=diff(ZZf);
ZMf = [];
for ik=1:kzz-1
  ZMf(ik,1)=ZZf(ik)+0.5*(ZZf(ik+1)-ZZf(ik));
end

for YR=YR1:YR2
  yr=YR;
%  fmatu=sprintf('%s%3.3i_SEGreenlSh_xsct_dayUTSZ_%i.mat',pthmat,expt,YR);
%  
  UTSZ=struct;
  UTSZ = SCT;
  UTSZ.Title = 'U norm, m/s, SE Greenl shelf sections, Z levels, monthly mean';
  UTSZ.ZZlevels = ZZf;
  UTSZ.ZM = ZMf;

  cc=0; 
 for iday=id1:id2
    dj1=datenum(yr,1,1);
    dnmb = dj1+iday-1; 
    DV   = datevec(dnmb);

    pthbin = sprintf('/nexsan/archive/ARCc0.08_011/data/%4.4i/',expt,yr);
    if expt==012,
      pthbin=sprintf('/nexsan/hycom/ARCc0.04_012/data/%4.4i/',yr);
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
    AA = sub_sct_interp2z(SCT,U,V,T,S,ZZf,ZZh,HH,dH,DX,DY);    
%keyboard
    UTSZ.Unrm(cc,:,:) = AA.NormalU;
    UTSZ.Temp(cc,:,:) = AA.Temp;
    UTSZ.Saln(cc,:,:) = AA.Salin;
    UTSZ.Time(cc)     = dnmb;	

    fprintf('++++++>  Processed 1 record %8.5f min\n\n',toc/60);

%    if mod(cc,10)==0 & s_mat>0
%      fprintf('Saving %s\n',fmatu);
%      save(fmatu,'UTSZ');
%    end	

  end  % dday
  
%  if s_mat>0
%    fprintf('Saving %s\n',fmatu);
%    save(fmatu,'UTSZ');
%  end	
  
end;   % years    
    
% Plot section
Dst=UTSZ.Dist_origin*1e-3;
S=squeeze(UTSZ.Saln);
T=squeeze(UTSZ.Temp);
U=squeeze(UTSZ.Unrm);
Hb=UTSZ.Hbottom;

S = sub_fill_bottom_nans(S);
T = sub_fill_bottom_nans(T);
% Filter to get rid of z-level disc.
WT=[0.3, 0.3, 0.25;...
    0.3, 0.3, 0.3;...
    0.3, 0.3, 0.3];
WT=WT./(sum(sum(WT)));
[nw,mw]=size(WT);
di=floor(nw/2);
dj=di;

[mm,nn]=size(S);
for it=1:6
dmm=S;
for ii=2:nn-1
  for jj=2:mm-1
    i1=ii-di;
    i2=ii+di;
    j1=jj-dj;
    j2=jj+dj;
    i1=max([1,i1]);
    i2=min([nn,i2]);
    j1=max([1,j1]);
    j2=min([mm,j2]);
    A=S(j1:j2,i1:i2);
    df=sum(sum(A.*WT));
    dmm(jj,ii)=df;
  end
end
S=dmm;
end


for it=1:6
dmm=T;
for ii=2:nn-1
  for jj=2:mm-1
    i1=ii-di;
    i2=ii+di;
    j1=jj-dj;
    j2=jj+dj;
    i1=max([1,i1]);
    i2=min([nn,i2]);
    j1=max([1,j1]);
    j2=min([mm,j2]);
    A=T(j1:j2,i1:i2);
    df=sum(sum(A.*WT));
    dmm(jj,ii)=df;
  end
end
T=dmm;
end


for it=1:6
dmm=U;
for ii=2:nn-1
  for jj=2:mm-1
    i1=ii-di;
    i2=ii+di;
    j1=jj-dj;
    j2=jj+dj;
    i1=max([1,i1]);
    i2=min([nn,i2]);
    j1=max([1,j1]);
    j2=min([mm,j2]);
    A=U(j1:j2,i1:i2);
    df=sum(sum(A.*WT));
    dmm(jj,ii)=df;
  end
end
U=dmm;
end

nfg=1;
i1=max(find(Hb>=-400));
i2=max(find(Hb>-2600));
xl1=Dst(i1);
xl2=Dst(i2);
yl1=-1000;
yl2=0;
stl=sprintf('0.04HYCOM-012, S, %s, %s',datestr(UTSZ.Time),UTSZ.Name); 
sub_plot_xsctS(ZMf,Dst,S,U,Hb,stl,nfg,xl1,xl2,yl1,yl2); % plot similar to MVP plot
DV=datevec(UTSZ.Time);
if s_fig>0
  fgnm=sprintf('%sarc004_%3.3i_vertS_%s_%i%2.2i%2.2i',...
       pthfig,expt,sctnm,DV(1:3));
  fprintf('Saving %s\n',fgnm);
  bottom_text(btx,'pwd',1);
  print('-dpng','-r150',fgnm);
end

nfg=2;
sub_plot_xsctT(ZMf,Dst,T,U,Hb,stl,nfg,xl1,xl2,yl1,yl2);
if s_fig>0
  fgnm=sprintf('%sarc004_%3.3i_vertT_%s_%i%2.2i%2.2i',...
       pthfig,expt,sctnm,DV(1:3));
  fprintf('Saving %s\n',fgnm);
  bottom_text(btx,'pwd',1);
  print('-dpng','-r150',fgnm);
end

  
  
  
