% Vertical section of T/S from GDEM4
% Monthly
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear


imo=1;  % 00 - annual

pfld  = 'temp';
%pfld  = 'salin';
fld0  = pfld;
%xname = 'BerNatl'; % Bering - N.Atl
%xname = 'BeaufNatl'; % Beauf. Shelf - South Norway
%xname = 'BeaufIcel';
xname = 'BafnFram'; % section around Greenland from Baffin to Fram str 
%xname = 'LaptSea'; % Laptev Sea

%s_mat  = 0; % overridden if  s_extr=0
s_fig = 0;

%pthmat  = '/nexsan/people/takis/lydia/HYCOM/data_mat/'; 
%pthmat2 = '/Net/ocean/ddmitry/HYCOM/ARCc/data_mat/';
pthmat2 = '/nexsan/people/takis/lydia/HYCOM/data_mat/';
pthmat  = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_mat/';
pthfig  = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/fig_GDEM4/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthin   = '/Net/data/GDEM4/';  % climatology data with no land

btx='vert_TS_xsectGDEM.m';

%fmat    = sprintf('%sADPTH_decadal.mat',pthmat); 

% Extract data for:
phi0 = 50;  % southmost latitude

ftopo = sprintf('%sdepth_ARCc0.08_07.nc',pthtopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(HH);

% Convert longitudes to -180 +180 rrange
% piece inserted from xsectAO.m

elon=LON; alat=LAT;

A=elon;
[my,nx]=size(elon);
for i=1:nx
for j=1:my
  long=elon(j,i);
  dmm=long;
  while (abs(dmm)>180)
    if dmm<0
      dmm=dmm+360;
    else
      dmm=dmm-360;
    end;
  end;  % while dmm
  elon(j,i)=dmm;
end;  % for j
end;  % for i
clear A;
% Grid cell spacing
[DX,DY]=sub_dx_dy(LON,LAT);


% Sections
switch(xname);
  case('BeaufIcel');
Xv=[  -85.3602600097656
      -60.5440979003906
      -50.6139831542969
      -53.5484313964844];

Yv=[  63.8842163085938
      64.8982925415039
      59.2981948852539
      53.0625038146973];
% 
% Use HYCOM 0.08 indices since
% GDEM has been interpolated onto ARCc0.08 grid
  IJs=[479, 1687;
       928, 1250;
       1046 848;
       927 539];
  xl1=1;
  xl2=6735;
  yl1=-4500;
  yl2=0;
%  cs1=30;
%  cs2=35;  
  cs1=3.36;  % log scale S=28.79
  cs2=3.56;  % log scale S=35.16
  ct1=-2;
  ct2=5;

 case('BafnFram');
  IJs=[  558         998
         497         644
         505         343
         638         140
         977         473
         960         622
        1041         882
         995        1044];
  xl1=1;
  xl2=10500;
  yl1=-5000;
  yl2=0;
  cs1=3.36;  % log scale S=28.79
  cs2=3.56;  % log scale S=35.16
  ct1=-2;
  ct2=5;

end

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

f_map=0;
if f_map>0
  figure(10); clf;
  hold on
  contour(HH,[0 0],'k');
  contour(HH,[-5000:1000:-100],'Color',[0.7 0.7 0.7]);
  plot(IIs,JJs,'b.-');
%  set(gca,'xlim',[100 600],...
%         'ylim',[200 1000]);
  axis('equal');
  title(sprintf('GDEMv4, Section %s',xname));
%  keyboard
end

fprintf('Section: %s, Flag fig: %i\n',xname,s_fig);


%-------------------------------
finT=sprintf('%sptgdemv4f%2.2i.nc4',pthin,imo);  % ptgdemv4f##.nc4
%finS=sprintf('%ssgdemv4f%2.2i.nc4',pthin,imo);

ZZ = -nc_varget(finT,'Depth');
%TT = nc_varget(finT,'Potential_Temperature');
%SS = nc_varget(finS,'salinity');
%LNg = nc_varget(finT,'Longitude');
%LT0 = nc_varget(finT,'Latitude');
%iy   = max(find(LT0<=phi0))-1;
%LTg  = LT0(iy+1:end);
%n    = length(LNg);
%m    = length(LTg);
%k    = length(ZZ);
%TT   = TT(:,iy+1:end,:);
%SS   = SS(:,iy+1:end,:);
   

fgdm=sprintf('%sGDEM_topo.mat',pthmat);
f_intrp_Hgdem=0; % =0 - upload interpoalted, =1- interpolat, very slow

if f_intrp_Hgdem>0
  Hgdem = sub_intrp_Hgdem(LAT,LON,HH,LNg,LTg);
% Quick nearest neighbor interpolation:
  fprintf('Saving GDEM topo %s\n',fgdm);
  save(fgdm,'Hgdem');
else
  fprintf('Loading GDEM topo %s\n',fgdm);
  load(fgdm); % interpoalted for 
end


% loading data
fmatT=sprintf('%sTave_GDEM_%2.2i.mat',pthmat2,imo)
fmatS=sprintf('%sSave_GDEM_%2.2i.mat',pthmat2,imo)
%fmatDP=sprintf('%sDPave_GDEM_00_full.mat',pthmat)
if strncmp(pfld,'temp',4)
  load(fmatT,'Tav');
else
  load(fmatS,'Sav');
end
%load(fmatDP,'DPav');


% 
% Find section
% Depth along transsect
for cntr=1:length(IJs)
		ii=IJs(cntr,1);
		jj=IJs(cntr,2);
		Ibtm(cntr)=min(find(ZZ<HH(jj,ii)));
end;

clear AA ZA
for cntr=1:length(IJs)
		ii=IJs(cntr,1);
		jj=IJs(cntr,2);

  if strncmp(pfld,'temp',4)
		  AA(:,cntr)=squeeze(Tav(:,jj,ii));
  else
		  AA(:,cntr)=squeeze(Sav(:,jj,ii));
  end
		ZA(:,cntr)=ZZ;

		AA(Ibtm(cntr):end,cntr)=nan;
%		ZA(Ibtm(cntr):end,cntr)=nan;

end;

Nintf=size(ZA,1);
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
				XL(1:Nintf,ii)=sum(dXX(1:ii));
		end
% Note correct distances: dXX is distance for "zigzag" sections
end   

stl=sprintf('GDEMv4, %s, Month=%i,  %s',...
												fld0,imo,xname);
nf=1;
xdr = 1;
f_layer=0;
[ma,na]=size(AA);
s0=[];
if strncmp(pfld,'salin',4)
% Convert to log scale for better depiction of high-value S
		A0=AA;

		s0=35.6;
    s0 = max([s0,max(max(A0))]);
		flg=1;
		AA=salin2exp(s0,A0,flg);
		c1=min(min(AA));
%		c2=max(max(AA));
    c2=5.6;
		cntrS=[34.8:0.05:35];  
		cntr=exp((35.6-cntrS).^(-1));
else
		c1=ct1;
		c2=ct2;
		cntr=[-1.:0.5:4];  
end

sub_plot_xsection3(nf,XL,ZA,AA,stl,fld0,f_layer,...
																			xl1,xl2,yl1,yl2,c1,c2,cntr,s0);
bottom_text(btx,'pwd',1);





