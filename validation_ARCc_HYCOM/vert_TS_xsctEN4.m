% Vertical section from EN4 gridded product
%

addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

pfld  = 'salin';
%pfld  = 'temp';
fld0  = pfld;
xname = 'BeaufIcel';
dnmb  = datenum(2015,7,1); % date to plot
DV1   = datevec(dnmb);


%pthmat  = '/nexsan/people/takis/lydia/HYCOM/data_mat/'; 
%pthmat2 = '/Net/ocean/ddmitry/HYCOM/ARCc/data_mat/';
pthdat  = '/Net/kronos/ddmitry/EN4_v4.2.1/2015_anls/';
pthmat  = '/Net/tholia/ddmitry/hycom/ARCc0.08/data_mat/';
pthmat2 = '/Net/ocean/ddmitry/HYCOM/ARCc0.08/112/data_mat/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthin   = '/Net/data/GDEM4/';  % climatology data with no land

btx='vert_TS_xsectEN4.m';

%
% Get coordinates of the section line saved in vert_TS_xsect008.m
fout=sprintf('%sarc008_xsect_%s.mat',pthmat2,xname);
fprintf('Loading %s\n',fout);
XYs=load(fout);

% Section:
switch(xname);
  case('BeaufIcel');
%
%  cs1=30;
%  cs2=35;  
  cs1=3.36;  % log scale S=28.79
  cs2=3.56;  % log scale S=35.16
  ct1=-2;
  ct2=5;
end


yr=DV1(1);
imo=DV1(2);
fin = sprintf('%sEN.4.2.1.f.analysis.g10.%4.4i%2.2i.nc',pthdat,DV1(1),DV1(2));
Z=nc_varget(fin,'depth');
LT=nc_varget(fin,'lat');
LN=nc_varget(fin,'lon');
if strncmp(pfld,'salin',4)
  F = squeeze(nc_varget(fin,'salinity'));
else
  F = squeeze(nc_varget(fin,'temperature')); % K
  F = F-273.15;
end


% Find indices of the section:
II=find(LN>180);
LN(II)=LN(II)-360;

Xv=XYs.Xl;
Yv=XYs.Yl;
nv=size(Xv,1);
cc=0;
for in=1:nv
  x0=Xv(in);
  y0=Yv(in);

  DI=abs(LN-x0);
  ii=find(DI==min(DI),1);
  DJ=abs(LT-y0);
  jj=find(DJ==min(DJ),1);


  if cc==0
    cc=cc+1;
    IJs(cc,1)=ii;
    IJs(cc,2)=jj;
  else  % check previous point
    di=abs(IJs(cc,1)-ii);
    dj=abs(IJs(cc,2)-jj);
    if di>0 | dj>0
      cc=cc+1;
      IJs(cc,1)=ii;
      IJs(cc,2)=jj;
    end
  end

end
nS=cc;

f_map=0;
if f_map==1
  figure(10); clf;
  a=squeeze(F(1,:,:));
  pcolor(a); shading flat;
  hold on;
  plot(IJs(:,1),IJs(:,2),'k.','Markersize',10);
end


%
% Distances:
for ii=1:nS
		i0=IJs(ii,1);
		j0=IJs(ii,2);
		Xl(ii)=LN(i0);
		Yl(ii)=LT(j0);
end

x1=Xl(1);
y1=Yl(1);
for ii=1:nS
		x2=Xl(ii);
		y2=Yl(ii);
		dXX(ii)=distance_spheric_coord(y1,x1,y2,x2)*1e-3;  % spheric distance, km
		x1=x2;
		y1=y2;
end
XL1=cumsum(dXX);

% Get section:
clear AA LTs LNs
for ii=1:nS
  i0=IJs(ii,1);
  j0=IJs(ii,2);
  AA(:,ii)=squeeze(F(:,j0,i0));
  LTs(ii)=LT(j0);
  LNs(ii)=LN(i0);
end

dll=diff(LNs);
I0=find(abs(dll)>100); % North Pole - no data


stl=sprintf('EN4v4.2.1 %s, %4.4i-%2.2i, %s',...
												fld0,yr,imo,xname);
nf=1;
xdr = 1;
f_layer=0;
[ma,na]=size(AA);
s0=[];
if strncmp(pfld,'salin',4)
% Convert to log scale for better depiction of high-value S
		A0=AA;

		s0=35.6;
		flg=1;
		AA=salin2exp(s0,A0,flg);
%		c1=min(min(AA));
%		c2=max(max(AA));
  c1=1.1316;  % to match arc008
  c2=5.56293;
		cntrS=[34.8:0.05:35];
		cntr=exp((35.6-cntrS).^(-1));
else
		c1=ct1;
		c2=ct2;
		cntr=[-1.:0.5:4];
end

% Missing data around the N Pole
%AA(:,I0)=nan;

xl1=1;
xl2=max(max(XL1));
yl1=-4500;
yl2=0;

dmm=ones(length(XL1),1);
ZZ=-Z*dmm';
dmm=ones(length(Z),1);
XL=dmm*XL1;

sub_plot_xsection3(nf,XL,ZZ,AA,stl,fld0,f_layer,...
																			xl1,xl2,yl1,yl2,c1,c2,cntr,s0);
% Small patch over pole region:
patch([XL1(I0-4) XL1(I0-4) XL1(I0+4) XL1(I0+4)],[-4500 0 0 -4500],'w');



bottom_text(btx,'pwd',1);





