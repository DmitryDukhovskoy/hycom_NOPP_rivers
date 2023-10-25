% Plot monthly climatology mean UV fields
% 1 month at a time
% with standard deviation ellipses at selected locations
% derived in mnthly_arc08_UV.m
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

s_fig  = 0;
imo = 7; % plot month=1,...,12 - 1 month at a time, 


plr =1;  % U from plr layer
rg = 9806;


regn = 'ARCc0.08';
expt = 110;
%pthfig  = '/nexsan/people/ddmitry/Net_ocean/hycom/ARCc0.08/110/fig_meanUV/';
pthfig  = '/Net/mars/ddmitry/hycom/ARCc0.08/110/fig_meanUV/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat =sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);

YRPLT = [1993:2016];
%YRPLT = [2013];
np = length(YRPLT);

fprintf('Month Climatology: %i\n',imo);

ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY] = sub_dx_dy(LON,LAT);

% Select locations where to plot StDev Ell.
% Arctic Ocean and N. Atlantic
% may need to adjust scales
% as speed is very different in 
% N. Atl and Arctic O. 
IJp =[631   454
   651   452
   679   537
   688   493
   790   543
   782   595
   857   597
   867   639
   913   662
   896   711
   940   713
   785   443
   657   390
   870   448
   826   526
   734   492];
npp = length(IJp);

usm=zeros(mm,nn);
vsm=zeros(mm,nn);
for ik=1:np
  iyr = YRPLT(ik);
  fmat = sprintf('%smnthUV_lr%2.2i_%i.mat',pthmat,plr,iyr);
  fprintf('Loading %s\n',fmat);
  load(fmat);
  
  U = meanUV(imo).U;
  V = meanUV(imo).V;
  usm = usm+U;
  vsm = vsm+V;

% Save data for Calculating StDev
  for ipp = 1:npp
    i0=IJp(ipp,1);
    j0=IJp(ipp,2);
    
    Up(ik,ipp)=meanUV(imo).U(j0,i0);
    Vp(ik,ipp)=meanUV(imo).V(j0,i0);
    
  end
  
  
end
U=usm/np;
V=vsm/np;

xlim1 = 620;
xlim2 = 980;
ylim1 = 370;
ylim2 = 770;

S = sqrt(U.^2+V.^2);
nf = 1;
stl = sprintf('ARCc0.08-%i, U, mo=%i, %i-%i',...
	      expt,imo,min(YRPLT),max(YRPLT));
c1=0;
c2=0.2;
f_cmp=3;
sub_plot_scalar_v2(S,nf,HH,xlim1,xlim2,...
		   ylim1,ylim2,LON,LAT,...
		   stl,'c1',c1,'c2',c2,'cmp',f_cmp);



scl=1.1;
cf=0.15;
beta=20;
lwd=1.2;
v_col=[0.3 0.3 0.3];
dii=5;
for ii=xlim1:dii:xlim2
  for jj=ylim1:dii:ylim2
    clear u v
    u = U(jj,ii)*100; % m/s -> cm/s
    v = V(jj,ii)*100;
    s = S(jj,ii);
    if isnan(u); continue; end;
    %    if res>0,
%    u=u/s;  % scale to speed, i.e. unit vectors
%    v=v/s;
%    scl=25;
    %    end

    x0=ii;
    y0=jj;

    x1=x0+u*scl;
    y1=y0+v*scl;
    draw_arrowF(x0,x1,y0,y1,cf,beta,v_col,lwd);
  end
end

% Calculate StDev ellipses
% using monthly data
fprintf('Plotting StDev ellipses ...\n');
um2cm=100;
%sclv=20;
sclv=scl;
for ipp = 1:npp
  clear uu vv
  i0=IJp(ipp,1);
  j0=IJp(ipp,2);
  if i0>=xlim1 & i0<=xlim2 ...
     j0>=ylim1 & j0<=ylim2
    uu=Up(:,ipp);
    vv=Vp(:,ipp);
    sub_std_ellipse(uu,vv,i0,j0,um2cm,sclv);
  end
  
end

f_legend=1;
if f_legend>0
  ulg=0.2*um2cm;
  vlg=0;
  ii=700;
  jj=650;
  cf=0.35;
  beta=20;
  col=[1 0 0];
  lwd=2.;
  draw_arrowF(ii,ii+ulg*sclv,jj,jj+vlg*sclv,cf,beta,col,lwd);
  ptx=text(ii,jj+10,sprintf('u=%3.2f m/s',ulg/um2cm));
  set(ptx,'Color',[1 0 0],'Fontsize',18);
end


%set(gca,'xlim',[280 1600],...
%	'ylim',[10 1950]);

txtb = 'plot_meanUVmonth_stdv_eastGr.m';
bottom_text(txtb,'pwd',1,'fontsize',8);

if s_fig>0
  fgnm = sprintf('%sarc08_110_meanUVmonth_std_Lr%2.2i_%i',...
		 pthfig,plr,iyr);
  fprintf('Saving %s\n',fgnm);
  print('-dpng','-r350',fgnm);
end




