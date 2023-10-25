% Plot T and U vectors, Lofoten basin
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

expt = '011';
TV = 17;
rg=9806;  % convert pressure to depth, m
huge=1e10;
omg = 7.2921159e-5; 

s_fig = 0;
zz0   = -100;  %  depth of calculation
%txtb = 'plot_deepU_ssh.m';
btxt = 'plot_Lofoten004.m';

pthmat  = '/Net/mars/ddmitry/hycom/ARCc0.08/data_mat/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/topo_grid/';
pthfig  = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/011/fig_Lofoten/';
%monmat = sprintf('%sEKE_meanUV_lev%i.mat',pthmat,Nlev);

ftopo = sprintf('%sdepth_ARCc0.04_17DD.nc',pthtopo); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

[DX,DY]=sub_dx_dy(LON,LAT);
Fcor = 2*omg*sind(LAT);
Lmsk = HH*0;
Lmsk(HH<zz0) = 1;

cnt = 0;
VRT = zeros(mm,nn);
Umn = zeros(mm,nn);
Vmn = zeros(mm,nn);

%
c1=0; 
c2=8;
%CMP = create_colormap5(200,0,8);
CMP = colormap_sclr2(200,c1,c2);
cmp=CMP.colormap;
cmp = smooth_colormap(cmp,18);

% Lofoten:
% 66-73 N latitude
% 10W - 18E long.
d1=distance_spheric_coord(66,-10,LAT,LON);
[j1,i1]=find(d1==min(min(d1)));
d1=distance_spheric_coord(73,-10,LAT,LON);
[j2,i2]=find(d1==min(min(d1)));
d1=distance_spheric_coord(73,18,LAT,LON);
[j3,i3]=find(d1==min(min(d1)));
d1=distance_spheric_coord(66,18,LAT,LON);
[j4,i4]=find(d1==min(min(d1)));
xl1=2060;
xl2=2460;
yl1=1160;
yl2=1500;
%figure('Visible','off'); clf;
figure(1); clf;

for yr=2005:2005
  fmat = sprintf('%smean_vort008_%i.mat',pthmat,yr);
  for iday=131:1:131
    tic;
    dnmb=datenum(yr,1,1)+iday-1;
%dnmb = datenum(2005,8,21); % date to plot
    DV =datevec(dnmb);
    yr = DV(1);
%iday = dnmb-datenum(yr,1,1)+1;
    pthbin = sprintf('/nexsan/hycom/ARCc0.04_%s/data/%4.4i/',expt,yr);

    fina = sprintf('%s%s_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
    finb = sprintf('%s%s_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);

    if ~exist(fina,'file') | ~exist(finb,'file');
      fprintf('Not found *a or *b: %s\n',fina);
      fprintf('                     %s\n',finb);
    %  continue;
    end

    %cnc=cnc+1;
    %TM(cnc,1)=dnmb;
    if ~exist('ZZ','var')
      [ZM,ZZ] = sub_zz_zm(fina,finb,HH);
      ZZ(isnan(ZZ))=100;
      ZM(isnan(ZM))=100;
      nlr = size(ZM,1);
      
      Iz=find(HH<zz0);
      for k=1:nlr
	zav(k) = nanmean(ZM(k,Iz));
      end
      dzav = abs(zav-zz0);
      Iz0 = find(dzav==min(dzav));
    end

    fprintf('Getting data expt %s: %4.4i_%2.2i_%2.2i: %s\n',expt,DV(1:3),fina);

    [F,n,m,nlr] = read_hycom(fina,finb,'u-vel.');
    %  F=F(:,jnc1:jnc2,inc1:inc2);
    F(F>huge)=0;
    U=squeeze(F(Iz0,:,:));

    [F,n,m,nlr] = read_hycom(fina,finb,'v-vel.');
    F(F>huge)=0;
    V=squeeze(F(Iz0,:,:));

    [F,n,m,nlr] = read_hycom(fina,finb,'temp');
    F(F>huge)=0;
    T=squeeze(F(Iz0,:,:));
    T(HH>0)=nan;
    
    clf;
    pcolor(T); shading flat;
    hold on;
%    contour(HH,[0 0],'k');
    contour(HH,[-5000:500:-200],'Color',[0.9 0.9 0.9]);
    caxis([c1 c2]);
    colormap(cmp);
    
    contour(LAT,[40:5:88],'Color',[0.8 0.8 0.8]);
    contour(LON,[-40:10:40],'Color',[0.8 0.8 0.8]);
    
    axis('equal');
    set(gca,'xlim',[xl1 xl2],...
	    'ylim',[yl1 yl2]);
    set(gca,'xtick',[],...
	    'ytick',[],...
	    'Color',[0 0 0]);
    ch = colorbar;

    scl=45;
    cf=0.3;
    beta=15;
    lwd=1.2;
    v_col=[0 0 0];
    dii=3;
    for ii=xl1:dii:xl2
      for jj=yl1:dii:yl2
	clear u v
	if HH(jj,ii)>-200, continue; end;
	u = U(jj,ii);
	v = V(jj,ii);
%	s = S(jj,ii);
	if isnan(u), continue; end;

%	u=u/s;
%	v=v/s;

	x0=ii;
	y0=jj;

	x1=x0+u*scl;
	y1=y0+v*scl;
	draw_arrowF(x0,x1,y0,y1,cf,beta,v_col,lwd);
      end
    end
%    colorbar;

    stl=sprintf('ARCc0.04-%s, T&U, depth=%5.1fm, %4.4i/%2.2i/%2.2i',expt,zz0,DV(1:3));
    title(stl);

    bottom_text(btxt,'pwd',1,'fontsize',10);
    cnt=cnt+1;
    if s_fig==1
      fgnm = sprintf('%sarc04_110_Lofoten_TU_%4.4i',pthfig,cnt);
      fprintf('Saving %s\n',fgnm);
      print('-dpng','-r200',fgnm);
    end

    fprintf('1 record: %6.4f min\n\n',toc/60);
    
  end
end
