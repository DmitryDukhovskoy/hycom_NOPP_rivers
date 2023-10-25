% Plot U at 100m 
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

expt = '011';
TV = '17DD';
rg=9806;  % convert pressure to depth, m
huge=1e10;
omg = 7.2921159e-5; 
dnmb = datenum(2006,1,10); % date to plot

sfig = 0;
zz0   = -100;  %  depth of calculation
%txtb = 'plot_deepU_ssh.m';
btx = 'plot_U_004.m';

pthmat  = '/Net/mars/ddmitry/hycom/ARCc0.04/data_mat/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/topo_grid/';
%pthfig  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/%s/fig_2D/',...
%		  expt);
pthfig = '/Net/mars/ddmitry/hycom/ARCc0.04/011/fig_vort/';

ftopo = sprintf('%sdepth_ARCc0.04_%s.nc',pthtopo,TV); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

[DX,DY]=sub_dx_dy(LON,LAT);
Fcor = 2*omg*sind(LAT);
Lmsk = HH*0;
Lmsk(HH<zz0) = 1;

DV =datevec(dnmb);
yr = DV(1);
iday = dnmb-datenum(yr,1,1)+1;

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

% find AWI grid boundary:
lt1=75;
lt2=82.5;
ln1=-20;
ln2=20;

d1=distance_spheric_coord(lt1,ln1,LAT,LON);
[j1,i1] = find(d1==min(min(d1)),1);
d1=distance_spheric_coord(lt2,ln1,LAT,LON);
[j2,i2] = find(d1==min(min(d1)),1);
d1=distance_spheric_coord(lt2,ln2,LAT,LON);
[j3,i3] = find(d1==min(min(d1)),1);
d1=distance_spheric_coord(lt1,ln2,LAT,LON);
[j4,i4] = find(d1==min(min(d1)),1);
jL1=min([j1,j2,j3,j4]);
jL2=max([j1,j2,j3,j4]);
%jL2=2200;
iL1=min([i1,i2,i3,i4]);
iL2=max([i1,i2,i3,i4]);
iL2=2290;

iV = [iL1,iL1,iL2,iL2];
jV = [jL1,jL2,jL2,jL1];
[II,JJ]=meshgrid([1:nn],[1:mm]);
dmm = inpolygon(II,JJ,iV,jV);
IN  = find(dmm==1);
IOUT= find(dmm==0);

HH(IOUT)=nan;

fprintf('Getting data expt %s: %4.4i_%2.2i_%2.2i: %s\n',expt,DV(1:3),fina);

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

[F,n,m,nlr] = read_hycom(fina,finb,'u-vel.','r_layer',Iz0);
%  F=F(:,jnc1:jnc2,inc1:inc2);
F(F>huge)=0;
UN=squeeze(F);
UN(IOUT)=nan;

[F,n,m,nlr] = read_hycom(fina,finb,'v-vel.','r_layer',Iz0);
F(F>huge)=0;
VN=squeeze(F);
VN(IOUT)=nan;
VN(HH>zz0)=nan;

CMP = create_colormap3(200,-0.3, 0.3);
cmpS=CMP.colormap;
for k=1:10
  cmpS(k,:)=[1 1 1];
%  cl2(k,:)=[1,1,1];
end
cmpS = smooth_colormap(cmpS,15);

S = sqrt(UN.^2+VN.^2);

%figure('Visible','off');
figure(1); clf;
pcolor(S); shading flat;
colormap(cmpS);
hold on
caxis([0 0.5]);
contour(HH,[0 0],'w');
axis('equal');
set(gca,'Color',[0 0 0],...
        'xlim',[iL1 iL2],...
        'ylim',[jL1 jL2]);
% Plot vectors
fprintf('Plotting vectors ...\n');
[mm,nn]=size(UN);

scl=14;
cf=0.3;
beta=18;
lwd=1.2;
v_col=[0 0 0];
dii=7;
for ii=iL1:dii:iL2
  for jj=jL1:dii:jL2
    clear u v
    u = UN(jj,ii);
    v = VN(jj,ii);
    s = S(jj,ii);
    if isnan(u), continue; end;
  %    if res>0,
    u=u/s;
    v=v/s;
%    scl=25;
  %    end

    x0=ii;
    y0=jj;

    x1=x0+u*scl;
    y1=y0+v*scl;
    draw_arrowF(x0,x1,y0,y1,cf,beta,v_col,lwd);
  end
end
colorbar;

stl=sprintf('|U|, HYCOM 0.04 %i/%2.2i/%2.2i',DV(1:3));
title(stl);


%txtb = 'plot_U_004.m';
bottom_text(btx,'pwd',1);
if sfig>0
  fgnm=sprintf('%sarc004_%s_U_%i%2.2i%2.2i',pthfig,expt,DV(1:3));
  fprintf('Saving fig %s\n',fgnm);
  print('-dpng','-r300',fgnm);
end


