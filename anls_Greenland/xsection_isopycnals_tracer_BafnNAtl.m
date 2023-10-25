% HYCOM Binary dayli mean files
% Plot vertical distribution of 
% tracer concentration (kg/m3) nat. log scale
% and
% plot mean vertical isopycnal layers 
% across Baffin-Labrador-nAtl-Nordic Seas
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

s_mat=0; % =0 - extract data, no save; =1 extract & save, =2 - load saved
s_fig=0;
rg=9806;  % convert pressure to depth, m

nTr = 1; % what tracer
fld0= 'tracer'; % tracer only
plr = 20; % highlight this interface

%figure(1); clf;
%set(gcf,'Visible','off');

%YPLT=[2004,299;2005,351;2006,295;2007,354;2008,286;2009,236];
% Create YRPLT array of 
% dates for averaging
% If plot snapshots - 1 field at a time in YRPLT
YRPLT=[];
cc=0;
iyr2=2010;
for iyr=2016:2016
  for idd=365:365
    if idd==1, idd=2; end;
    cc=cc+1;
    YRPLT(cc,1)=iyr;
    YRPLT(cc,2)=idd;
  end
end

%YRPLT=[2008,50];
np=size(YRPLT,1);


regn = 'ARCc0.08';
expt = 110;
pthfig  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/fig_sections/',expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);

fprintf('Plotting xsection for %s-%3.3i, Tracer %2.2i\n\n',regn,expt,nTr);

ftopo = sprintf('%s/depth_%s_09.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

% Specify segments of the x-section
IJs=[    558         998
         497         644
         505         343
         638         140
         977         473
         960         622
        1041         882
         995        1044];

IIa=IJs(:,1);
JJa=IJs(:,2);

     
nij=size(IJs,1);
IIs=[];
JJs=[];
for ii=1:nij-1
  i1=IJs(ii,1);
  i2=IJs(ii+1,1);
  j1=IJs(ii,2);
  j2=IJs(ii+1,2);
  [I,J]=sub_xsct_indx(i1,j1,i2,j2);
  if size(I,1)==1;
    I=I';
    J=J';
  end
  
  IIs=[IIs;I];
  JJs=[JJs;J];
end;

IJs=[IIs,JJs];

nS=length(IIs);
clear Xl Yl
for ii=1:nS
  i0=IJs(ii,1);
  j0=IJs(ii,2);
  Xl(ii,1)=LON(j0,i0);
  Yl(ii,1)=LAT(j0,i0);
end;
INDs=sub2ind(size(HH),JJs,IIs);

% 
% Segment end points:
for ii=1:length(IIa)
  i0=IIa(ii);
  j0=JJa(ii);
  XXa(ii,1)=LON(j0,i0);
  YYa(ii,1)=LAT(j0,i0);
end

f_map=0;
if f_map>0
  figure(10); clf;
  hold on
  contour(HH,[0 0],'Color',[0.6 0.6 0.6]);
  contour(HH,[-5000:1000:-100],'Color',[0.8 0.8 0.8]);
  plot(IIs,JJs,'b.-');
%  set(gca,'xlim',[100 600],...
%   'ylim',[200 1000]);
  axis('equal');
  set(gca,'xlim',[300 1200],...
          'ylim',[100 1100]);

  for ii=1:length(IIa)
    i0=IIa(ii);
    j0=JJa(ii);

    x0=XXa(ii);
    y0=YYa(ii);
    plot(i0,j0,'r.','Markersize',14);
    lts=sprintf('%7.3fW %7.3fN',x0,y0);
    text(i0+5,j0,lts,'Fontsize',14);
  end

end

%keyboard

% Plot fields:
cnc=0;  % for averaging, uncomment and comment in the loop
ip1=1;
for ip=ip1:np
%  cnc=0;   % do not average
  yr=YRPLT(ip,1);
  iday=YRPLT(ip,2);

  pthbin = sprintf('/nexsan/archive/%s_%3.3i/data/%i/',regn,expt,yr)
  fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
  finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);


  if ~exist(fina,'file');
    fprintf('Not found: %s\n\n',fina);
    continue;
  end
  
  cnc=cnc+1;
  dnmb=datenum(yr,1,1)+iday-1;
  DV=datevec(dnmb);

  fprintf('%4.4i_%2.2i_%2.2i: %s\n',DV(1:3),fina);

  tic;
%  [F,n,m,l] = read_hycom(fina,finb,fld0);
  [F,n,m,l] = read_hycom(fina,finb,'tracer','r_tracer',nTr);
  toc;
  F(F>1e6)=nan;

  Tr=squeeze(F(:,INDs));
  Tr(Tr<=1e-12)=nan;
%  lTr=log(Tr);
% Prepare for plotting - add extra bogus layer
% at the bottom, Matlab won't show it
  [a1,a2]=size(Tr);
%  Tr(end+1,1:a2)=Tr(end,:);
  
  
% Prepare vertical thicknesses for plotting  
  fld='thknss';
  [F,n,m,l] = read_hycom(fina,finb,fld);
  F(F>1e10)=nan;
  F(F<0.1)=nan;
  F=F/rg;
  Dsec=squeeze(F(:,INDs));
  Dsec(Dsec==0)=nan;
% Create Depth array of interface depths:
% Note these are BOTTOM interfaces 
% So Layer 1 is between interfaces 0m and ZZ(1)
  clear ZZb
  Dsec(isnan(Dsec))=0;
  ZZb(1,:)=-Dsec(1,:);
  for kk=2:l
    ZZb(kk,:)=ZZb(kk-1,:)-Dsec(kk,:);
  end
% For plotting need to have #of layers +1 interfaces
% add surface layer, otherwise all values
% will be shifted 1 layer down
  [nl,npb]=size(ZZb);
  ZZ=zeros(nl+1,npb);
  ZZ(2:nl+1,:)=ZZb;
  if cnc==1
    ZZav=ZZ*0;
    DsecAv=Dsec*0;
    TrZ=Dsec*0;
  end
  ZZav=ZZav+ZZ;
%   ZZav=ZZ;
  TrZ=TrZ+Tr.*Dsec;  % Field * layer thickness - for averaging
  DsecAv=DsecAv+Dsec;
%  if ~exist('Dav','var')
%    Dav=Dlr*0;
%  end
%  Dav=Dav+Dlr;
%
% Depths of middle of the cells:
  ZM(1,:)=0.5*ZZ(1,:);
  for kk=1:l
    ZM(kk,:)=0.5*(ZZ(kk+1,:)+ZZ(kk,:));
  end

 
  Nintf=size(ZZ,1);
  if ~exist('XL','var')
    nn=length(Xl);
    x1=Xl(1);
    y1=Yl(1);
    for ii=1:nn
      x2=Xl(ii);
      y2=Yl(ii);
      dXX(ii)=distance_spheric_coord(y1,x1,y2,x2)*1e-3;  % spheric distance, km
      x1=x2;
      y1=y2;
      XL(1:Nintf,ii)=sum(dXX);
    end
% Correct distances: dXX is distance for "zigzag" sections
% that follow the straight lines, thus the total
% distancs is longer than the straight line distances
% need to adjust it:
%    xtot=7922; % that what it should be, km, from Myers
%    xmx=max(max(XL));
%    crct=xtot/xmx;
%    XL=XL*crct;
  end   


end;

%keyboard

DsecAv(DsecAv==0)=nan;
Trav=TrZ./DsecAv;  % mean field (salinity or tracer conc.)
ZZav=ZZav/cnc;
% Add extra layer at ZZ=0 for plotting:
Trav=[Trav(1,:);Trav(:,:)];
%  Dav=Dav/cnc;

lTr=Trav;
if fld0=='tracer'
%  Trav(Trav<=1e-3)=1e-20;
%  Trav(isnan(Trav))=1e-20;
  lTr=log(Trav);
end  


%%  lt0=mean(LAT(j1,i1:i2));
ln1=LON(j1,i1);
ln2=LON(j2,i1);

stl=sprintf('Tr %i %4.4i/%2.2i/%2.2i, BLayer=%i,BaffinB-NAtl-Nordic Seas',...
	    nTr,DV(1:3),plr);
nf=1;
f_layer=0; % plot sigma-layers or not
sub_plot_xsection(nf,XL,ZZav,lTr,stl,fld0,f_layer);
btx = 'xsection_isopycnals_tracer_BafnNatl.m';
bottom_text(btx,'pwd',1);

%  ffg=sprintf('%stracers_conc_%4.4i_%2.2i_%2.2i-layer%2.2i',pthfig,DV(1:3),plr);
ffg=sprintf('%sxsct_BaffNAtlNord-Tr%2.2i-%4.4i%2.2i%2.2i',...
	    pthfig,nTr,DV(1:3));

%  keyboard

if s_fig>0
  fprintf('Saving %s\n\n',ffg);
  print('-dpng','-r200',ffg);
end

% Plot the cross-section
LMSK = HH*0;
LMSK(HH<0)=1;

cbw = [0 0 0; 1 1 1];

figure(2); clf;
hold on
pcolor(LMSK); shading flat;
colormap(cbw);
freezeColors;
contour(HH,[0 0],'k');
contour(HH,[-5000:1000:-100],'Color',[0.7 0.7 0.7]);
plot(IIs,JJs,'b.');
x=XL(1,:);
for kk=1000:1000:max(x)
  d=abs(x-kk);
  i0=find(d==min(d));
  plot(IIs(i0),JJs(i0),'ro','linewidth',2);
end
[mm,nn]=size(HH);
axis('equal');
set(gca,'xlim',[200 1200],...
	'ylim',[1 1200]); 
set(gca,'xtick',[],'ytick',[]);
btx = 'xsection_isopycnals_tracer_BafnNatl.m';
bottom_text(btx,'pwd',1);

ffg=sprintf('%ssgm_BaffNAtlNord-map',pthfig);

if s_fig>0
  fprintf('Saving %s\n\n',ffg);
  print('-dpng','-r200',ffg);
end



