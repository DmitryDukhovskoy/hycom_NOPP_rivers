% Extract T/S at Arctic Shelves
% Daily by months with time step dday
% 0.04 HYCOM-CICE
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

YR=2017;
mo1=6;
mo2=6;
dmm=1; % month step
dday=7; 


s_mat=1; 

rg = 9806;
rhow=1027; 
hgg=1e20;

regn = 'ARCc0.04';
%expt = 011;
%expt = 012; % Greenland runoff
expt = 022; % GOFS3.5

pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';
pthmat = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/%3.3i/data_shelf/',expt);
btx='shelf_xsct004.m';

fprintf('%s-%3.3i Daily Shelf SE Section T,S Zlevels, %i/%i-%i/%i\n\n',...
	regn,expt,YR,mo1,YR,mo2);

ftopo = sprintf('%s/depth_%s_17DD.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);

i=0;
% Sections:
i=i+1;
SCT(i).Name='Beauf';
SCT(i).IJ=[953 3354; 989 3303];

i=i+1;
SCT(i).Name='EastSib';
SCT(i).IJ=[1716, 3289;...
           1957, 3649];

i=i+1;
SCT(i).Name='Laptev';
SCT(i).IJ=[2343, 3119;...
           2458, 3320];

i=i+1;
SCT(i).Name='Kara';
SCT(i).IJ=[2321, 2484; ...
           2958, 2484];

nsct=length(SCT);


for ip=1:nsct
  fprintf('Section: %s\n',SCT(ip).Name);
  IJs=SCT(ip).IJ;
  [IIs,JJs]=sub_xsct_indx(IJs(1,1),IJs(1,2),IJs(2,1),IJs(2,2));
  SCT(ip).I=IIs;
  SCT(ip).J=JJs;
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
  figure(fn); clf;
  hold on;
  contour(HH,[0 0],'k');
  contour(HH,[-5000:1000:-10],'Color',[0.7 0.7 0.7]);

  for ip=1:nsct
%    IJ=SCT(ip).IJ;
%    plot([IJ(1,1) IJ(2,1)],[IJ(1,2) IJ(2,2)],...
%	 'Linewidth',2.5,'Color',[1. 0.6 0]);
    IIs=SCT(ip).I;
    JJs=SCT(ip).J;
    plot(IIs,JJs,'-',...
	 'Linewidth',2.5,'Color',[1. 0.4 0]);
  end
  axis('equal');
  set(gca,'xlim',[800 3200],...
	  'ylim',[1500 3800]);
  bottom_text(btx,'pwd',1);
end

% Interpolate onto fixed z-grid
% Specify fixed ZZ (interf) and ZM levels:
ZZf = [(0:-1:-50)';(-52:-2:-100)';(-105:-5:-500)';...
       (-510:-10:-1200)';(-1500:-500:-5000)'];
kzz = length(ZZf);

dZf=diff(ZZf);
ZMf = [];
for ik=1:kzz-1
  ZMf(ik,1)=ZZf(ik)+0.5*(ZZf(ik+1)-ZZf(ik));
end


for MO=mo1:dmm:mo2;
  dnmb1=datenum(YR,MO,1);
  dv=datevec(dnmb1+32);
  dnmb2=datenum(dv(1),dv(2),1)-1;

  fmatu=sprintf('%s%3.3i_ArctShelf_xsct_dayTSZ_%4.4i%2.2i.mat',pthmat,expt,YR,MO);

  TSZ = SCT;
  for ip=1:nsct
    TSZ(ip).Title = 'T&S Arctic shelf sections, Z levels, daily';
    TSZ(ip).ZZlevels = ZZf;
    TSZ(ip).ZM = ZMf;
  end

  cc=0;
  for dnmb=dnmb1:dday:dnmb2
    DV=datevec(dnmb);
    yr=DV(1);
    dj1=datenum(yr,1,1);
    iday=dnmb-dj1+1;

%    pthbin = sprintf('/nexsan/archive/ARCc0.08_011/data/%4.4i/',expt,yr);
    pthbin = sprintf('/nexsan/people/ddmitry/hycom/ARCc0.04_022/data/%i_mean/',yr);
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

    [ZMh,ZZh] = sub_zz_zm(fina,finb,HH);
    dH=abs(diff(ZZh,1));
    AA = sub_sct_TSinterp2z(SCT,T,S,ZZf,ZZh,HH,dH,DX,DY);    
  %keyboard

    for ip=1:nsct
      TSZ(ip).Temp(cc,:,:) = AA(ip).Temp;
      TSZ(ip).Saln(cc,:,:) = AA(ip).Salin;
      TSZ(ip).Time(cc)     = dnmb;	
    end

    fprintf('++++++>  Processed 1 record %8.5f min\n\n',toc/60);

  end;   % days

  if s_mat>0
    fprintf('Saving %s\n',fmatu);
    save(fmatu,'TSZ');
  end	

end
