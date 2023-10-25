% I. Yashayaev data 
% surface u from drifters
% SPNA
addpath /Net/Movies0/ddmitry/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /home/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers;
%addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
%addpath /Net/Movies0/ddmitry/MyMatlab/seawater
startup

close all
clear

pthdat = '/nexsan/people/ddmitry/Net_ocean/drifters_yash/';
pthtopo= '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/110/data_GG_prt/';


ftopo = sprintf('%sdepth_ARCc0.08_11.nc',pthtopo); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);

load('SPG_noNorth_indx.mat');
IGR(end+1,:)=IGR(1,:);

% Get indices of segments
clear SCT
nsct=length(IGR)-1;
SCT.IGR=IGR;
for ii=1:nsct
  SCT(ii).IJ = [IGR(ii,:); IGR(ii+1,:)];
  IJs=SCT(ii).IJ;
  [IIs,JJs]=sub_xsct_indx(IJs(1,1),IJs(1,2),IJs(2,1),IJs(2,2));
  SCT(ii).I=IIs;
  SCT(ii).J=JJs;
  nsg=length(IIs);
  clear XX YY Hb
  for ij=1:nsg
    i0=IIs(ij);
    j0=JJs(ij);
    XX(ij)=LON(j0,i0);
    YY(ij)=LAT(j0,i0);
    Hb(ij)=HH(j0,i0);
  end
  SCT(ii).long=XX;
  SCT(ii).latd=YY;
  SCT(ii).Hb=Hb;
  iocn = length(find(Hb<0));
  SCT(ii).Nocn = iocn; % # of ocean points in the segment
end

btx = 'anls_drifters.m';
f_map=0;
if f_map==1
  figure(10); clf;
  hold on;

  lcmp=[0 0 0; 1 1 1];
  Lmsk=HH*0;
  Lmsk(HH<0)=1;
  pcolor(Lmsk); shading flat;
  colormap(lcmp);

  contour(HH,[-5000:500:0],'Color',[0.8 0.8 0.8]);
  caxis([0 1]);

  for isg=1:nsct
    II = SCT(isg).I;
    JJ = SCT(isg).J;

    iocn = SCT(isg).Nocn;
    if iocn<=0
      clr = [0.7 0.7 0.7];
    else
      clr = [1 0.4 0];
    end

    II = SCT(isg).I;
    JJ = SCT(isg).J;

    plot(II,JJ,'.-','Color',clr);
  end

  axis('equal');
  set(gca,'xlim',[300 1150],...
            'ylim',[50 800]);

  title('SPNA region');
  bottom_text(btx,'pwd',1);

  fsct = sprintf('%sSCT_SPNA.mat',pthmat);
  fprintf('Saving Sections SPNA %s\n',fsct);
  save(fsct,'SCT');

keyboard
end

% 
% Load drifters data
file_name = 'Mean, Year range 0-9999, Seasonal range 1, Region 34.csv';
fin=sprintf('%s%s',pthdat,file_name);
%NM = csvread(fin,0,0,[0 0 1 8]);
A = csvread(fin,1,0);
X=A(:,1);
Y=A(:,2);
U=A(:,3);
V=A(:,4);

y1=min(Y);
y2=max(Y);
x1=min(X);
x2=max(X);
dx=0.5;
dy=dx;

XX=[x1:dx:x2];
YY=[y1:dy:y2];

nx=length(XX);
ny=length(YY);

for ii=1:nx
  for jj=1:ny
    x0=XX(ii);
    y0=YY(jj);
    II=find(X==x0 & Y==y0);

    if isempty(II)
      UU(jj,ii)=nan;
      VV(jj,ii)=nan;
    else
      UU(jj,ii)=U(II);
      VV(jj,ii)=V(II);
    end
  end
end
s=sqrt(UU.*UU+VV.*VV);











