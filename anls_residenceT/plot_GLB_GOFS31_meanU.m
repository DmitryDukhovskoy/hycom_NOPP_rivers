% from Global HYCOM reanalysis GOFS3.1 (expt_53.x)
%
% SPG region
%
% Plot monthly climatology mean UV fields
% Can do 1 month at a time:
% im1=im2 for all years (monthly climatology)
% or average over all months
% im1:im2 if im1/=im2
%
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

s_mat  = 0; % =0 - load saved,  =1 - save
im1 = 1; % im1=im2 plot month=1,...,12 - 1 month at a time, 
im2 = 1; % otherwise - average over N years and all months


rg=9806;  % convert pressure to depth, m
Sref=35; % N.Atl. is too saline

pthARC = '/Net/mars/ddmitry/hycom/GLBb2ARCc0.08/mnth_mean/';
pthmat = '/Net/mars/ddmitry/hycom/GLBb0.08/data_mat/';
pthfig = '/Net/mars/ddmitry/hycom/GLBb0.08/fig_ssh/';
%ftopo = sprintf('%sGLBb_T07_subset_Natl.mat',pthmat);
%pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
%pthglb = '/nexsan/GLBb0.08/GLBb0.08_191/topo/';

%fmat = sprintf('%sGLB_GOFS31_annualSSH_SPG.mat',pthmat);
fout = sprintf('%sGLB_GOFS31_meanUV50m_SPG_%2.2i_%2.2i-1994-2015.mat',pthmat,im1,im2);


i1=1210;
i2=2445;
ni=i2-i1;
j1=2100;
j2=3100;
nj=j2-j1;


% Select locations where to plot StDev Ell.
% same as in plot_meanUmonth_stdv_SPG.m
XYp=[-46.702    59.654
-50.264    61.580
-57.103    63.191
-61.172    60.222
-57.935    56.859
-51.223    53.273
-46.300    49.294
-38.442    50.831
-31.045    50.947
-26.798    51.458
-19.247    53.156
-16.497    55.404
-13.084    57.275
-11.019    60.021
-13.922    60.863
-18.682    62.319
-28.205    59.814
-32.793    58.200
-31.203    60.973
-30.567    64.702
-35.556    63.783
-40.920    60.812
-43.283    59.125
-52.071    59.994
-57.659    60.958
-56.247    58.207
-53.197    57.385
-49.941    56.863
-46.577    54.212
-42.131    56.609
-35.183    60.120
-34.140    62.341
-33.480    55.007
-39.387    54.382
-29.374    55.014
-23.323    60.012
-18.490    56.413];
npp = length(XYp);

fnm='/nexsan/hycom/GLBv0.08_53X/data/1994/hycom_GLBv0.08_530_1994010112_t000.nc';
lon=nc_varget(fnm,'lon',[i1],[ni]);
lat=nc_varget(fnm,'lat',[j1],[nj]);
[LON,LAT]=meshgrid(lon,lat);
[IC,JC] = sub_find_indxHYCOM(LON,LAT,XYp(:,1),XYp(:,2));
IJp=[IC,JC];


if s_mat==1
  ikk=0;
  usm=zeros(nj,ni);
  vsm=zeros(nj,ni);
  for iyr=1994:2015
    pthbin=sprintf('/nexsan/hycom/GLBv0.08_53X/data/%i/',iyr);
  %  fmat = sprintf('%sGLB_GOFS31_annualSSH_SPG_%4.4i.mat',pthmat,iyr);


    for imo=im1:im2
      for iday=1:7:31
	if iyr<1999
	 expt=530;
	elseif iyr==1999 & imo<4
	  expt=530;
	elseif iyr==1999 & imo>=4
	  expt=531;
	elseif iyr==2000
	  expt=531;
	elseif iyr>=2001 & iyr<2003
	  expt=532;
	elseif iyr==2003 & imo<=6
	  expt=532;
	elseif iyr==2003 & imo>6
	  expt=533;
	elseif iyr>2003 & iyr<2005
	  expt=533;
	elseif iyr==2005 & imo<=6
	  expt=533;
	elseif iyr==2005 & imo>6
	  expt=534;
	elseif iyr>2005 & iyr<2007
	  expt=534;
	elseif iyr==2007 & imo<=6
	  expt=534;
	elseif iyr>2007 & iyr<2009
	  expt=535;
	elseif iyr==2009 & imo<=6
	  expt=535;
	elseif iyr==2009 & imo>6
	  expt=536;
	elseif iyr>2009 & iyr<2011
	  expt=536;
	elseif iyr==2011 & imo<=6
	  expt=536;
	elseif iyr==2011 & imo>6
	  expt=537;
	elseif iyr>2011 & iyr<2013
	  expt=537;
	elseif iyr==2013 & imo<=6
	  expt=537;
	elseif iyr==2013 & imo>6
	  expt=538;
	elseif iyr==2014
	  expt=538;
	elseif iyr==2015
	  expt=539;
	end

	ihr=0;
	fnm=sprintf('%shycom_GLBv0.08_%3.3i_%4.4i%2.2i%2.2i12_t%3.3i.nc',...
		    pthbin,expt,iyr,imo,iday,ihr);

	if ~exist(fnm,'file');
	  while ihr<21
	    ihr=ihr+3;
	    fnm=sprintf('%shycom_GLBv0.08_%3.3i_%4.4i%2.2i%2.2i12_t%3.3i.nc',...
		    pthbin,expt,iyr,imo,iday,ihr);
	     if exist(fnm,'file'); break; end;
	  end
	end

	if ~exist(fnm,'file');
	  fprintf('Does not exist: %s\n',fnm);
	  continue;
	end

	fprintf('Reading %s\n',fnm);

	tic;
	ikk=ikk+1;
	
	nz=15;
	if ~exist('DZ','var')
          zz=nc_varget(fnm,'depth');
	  dz0=abs(diff(zz));
	  dz=dz0(1:nz-1);
	  zbt=sum(dz);
	  DZ=zeros(nz-1,nj,ni);
	  for iz=1:nz-1
	    DZ(iz,:,:)=dz(iz);
	  end
% Sea/land mask
	  HH=squeeze(nc_varget(fnm,'surf_el',[0 j1 i1],[1 nj ni]));
	  Hmsk=HH*0;
	  IW=find(~isnan(HH));
	  Hmsk(IW)=1;
	  IL=find(isnan(HH));
	  Hmsk(IL)=0;
	end
	
	dmm=squeeze(nc_varget(fnm,'water_u',[0 0 j1 i1],[1 nz-1 nj ni]));
% Depth average
        U=squeeze(nansum(dmm.*DZ,1))/zbt;
	U(IL)=nan;

	dmm=squeeze(nc_varget(fnm,'water_v',[0 0 j1 i1],[1 nz-1 nj ni]));
% Depth average
        V=squeeze(nansum(dmm.*DZ,1))/zbt;
	V(IL)=nan;
	
	usm = usm+U;
	vsm = vsm+V;
	
%S=sqrt(U.^2+V.^2);
    % Save data for Calculating StDev
      for ipp = 1:npp
	i0=IJp(ipp,1);
	j0=IJp(ipp,2);

	Up(ikk,ipp)=U(j0,i0);
	Vp(ikk,ipp)=V(j0,i0);
      end

      s=sqrt(usm.^2+vsm.^2);
      fprintf('Max mean |U|=%7.2f m/s\n',max(max(s))./ikk);
      fprintf('1 record processing: %5.3f min\n\n',toc/60);
      
      end   % iday
    end     % imo
  end;      % year
  
  U=usm/ikk;
  V=vsm/ikk;

  UMN.U=U;
  UMN.V=V;
  UMN.yr1=1994;
  UMN.yr2=2015;
  UMN.tser_IJ=IJp;
  UMN.Utser=Up;
  UMN.Vtser=Vp;
  UMN.Hmsk=Hmsk;
  fprintf('Saving %s\n',fout);
  save(fout,'UMN');

else
  fprintf('Loading %s\n',fout);
  load(fout);
  
end

Hmsk=UMN.Hmsk;
U=UMN.U;
V=UMN.V;
Up=UMN.Utser;
Vp=UMN.Vtser;

xlim1 = 180;
xlim2 = 970;
ylim1 = 95;
ylim2 = 590;

S = sqrt(U.^2+V.^2);
nf = 1;
stl = sprintf('HYCOM GLB0.08 GOFS3.1, U, mo=%i/%i, %i-%i',...
	      im1,im2,1994,2015);

figure(1); clf;
%hmsk=HH;
%hmsk(HH<0)=nan;
pcolor(Hmsk); shading flat;
hold on;
%colormap([0.1 0.1 0.1]);
colormap([0.8 0.8 0.8; 1 1 1]);
%freezeColors;

%contour(HH,[-8000:500:-10],...
%	'Color',[0.6 0.6 0.6],...
%	'Linewidth',1.2);
axis('equal');
set(gca,'xlim',[xlim1 xlim2],...
	'ylim',[ylim1 ylim2]);
set(gca,'xtick',[],'ytick',[]);
%clr=[0.9 0.9 0.9];
%plot_gridlines(45,10,1,clr,LON,LAT);
title(stl,'Fontsize',12,'Interpreter','none');



scl=1.1;
cf=0.35;
beta=20;
lwd=1.6;
v_col=[0.3 0.3 0.3];
dii=10;
for ii=xlim1:dii:xlim2
  for jj=ylim1:dii:ylim2
    clear u v
    u = U(jj,ii)*100; % m/s -> cm/s
    v = V(jj,ii)*100;
    s = sqrt(u*u+v*v);
    if isnan(u); continue; end;
    
    u=u/s;  % scale to speed, i.e. unit vectors
    v=v/s;
    if s>=1 & s<5
      scl=4;
      v_col=[0.5 0.5 0.5];
      cf=0.8;
      beta=30;
    elseif s>=5 & s<10
      scl=10;
      v_col=[0. 0.6 1];
      cf=0.4;
    elseif s>=10 & s<15
      scl=21;
      v_col=[0.9 0.5 0.];
      cf=0.3;
    elseif s>=15
      scl=25;
      v_col=[0.5 0. 0.];
      cf=0.3;
    end

    x0=ii;
    y0=jj;

    x1=x0+u*scl;
    y1=y0+v*scl;
    draw_arrowF(x0,x1,y0,y1,cf,beta,v_col,lwd);
    
%    keyboard
    
  end
end

x0=405;
y0=580;
dyy=15;
u=1; v=0;
% plot scale vectors:
%
% 1-5 cm
scl=4;
v_col=[0.5 0.5 0.5];
cf=0.8;
beta=30;
x1=x0+u*scl;
y1=y0+v*scl;
draw_arrowF(x0,x1,y0,y1,cf,beta,v_col,lwd);
text(x1+20,y1,'1-5 cm/s','Fontsize',14,'Color',v_col);

beta=20;
% 5-10 cm
scl=10;
v_col=[0. 0.6 1];
cf=0.4;
y0=y0-dyy;
x1=x0+u*scl;
y1=y0+v*scl;
draw_arrowF(x0,x1,y0,y1,cf,beta,v_col,lwd);
text(x1+20,y1,'5-10 cm/s','Fontsize',14,'Color',v_col);

%10-15
scl=21;
v_col=[0.9 0.5 0.];
%v_col=[0. 0. 0.6];
cf=0.3;
y0=y0-dyy;
x1=x0+u*scl;
y1=y0+v*scl;
draw_arrowF(x0,x1,y0,y1,cf,beta,v_col,lwd);
text(x1+10,y1,'10-15 cm/s','Fontsize',14,'Color',v_col);

%>15
scl=25;
v_col=[0.5 0. 0.];
cf=0.3;
y0=y0-dyy;
x1=x0+u*scl;
y1=y0+v*scl;
draw_arrowF(x0,x1,y0,y1,cf,beta,v_col,lwd);
text(x1+10,y1,'15 cm/s','Fontsize',14,'Color',v_col);




% Calculate StDev ellipses
% using monthly data
fprintf('Plotting StDev ellipses ...\n');
%um2cm=100;
%sclv=20;
sclv=1.5;
lwd=1.8;
v_col=[0 0.2 0.5];
for ipp = 1:npp
  clear u v
  i0=IJp(ipp,1);
  j0=IJp(ipp,2);
  if i0>=xlim1 & i0<=xlim2 & ...
     j0>=ylim1 & j0<=ylim2
    u=Up(:,ipp)*100;
    v=Vp(:,ipp)*100;
    
    sub_std_ellipse_v2(u,v,i0,j0,sclv,v_col,cf,lwd,beta);
  end
  
end

f_legend=1;
if f_legend>0
  ulg=30;
  vlg=0;
  ii=x0;
  jj=y0-dyy;
  col=v_col;
  [Xh,Yh]=get_arrowF(ii,ii+ulg*sclv,jj,jj+vlg*sclv,cf,beta,col,lwd);
  ptx=text(Xh(1)+10,Yh(1),sprintf('<U> %2.2i cm/s',ulg));
  set(ptx,'Color',col,'Fontsize',14);

  ulg=15;
  vlg=0;
  ii=x0;
  jj=jj-dyy;
  col=v_col;
  [Xh,Yh]=get_arrowF(ii,ii+ulg*sclv,jj,jj+vlg*sclv,cf,beta,col,lwd);
  ptx=text(Xh(1)+10,Yh(1),sprintf('<U> %2.2i cm/s',ulg));
  set(ptx,'Color',col,'Fontsize',14);
end


%set(gca,'xlim',[280 1600],...
%	'ylim',[10 1950]);

txtb = 'plot_GLB_GOFS31_meanU.m';
bottom_text(txtb,'pwd',1,'fontsize',7);





