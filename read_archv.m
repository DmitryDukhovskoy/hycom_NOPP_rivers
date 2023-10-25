% Read instant output archive files:
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

expt='100';
yr=1993;
iday=54;
hr=0;

j0=2134;
i0=896;


pthbin  = '/Net/kronos/ddmitry/hycom/arc0.08/output/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';

ftopo = sprintf('%sdepth_ARCc0.08_09.nc',pthtopo); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');

fina = sprintf('%s%s_archv.%4.4i_%3.3i_%2.2i.a',pthbin,expt,yr,iday,hr);
finb = sprintf('%s%s_archv.%4.4i_%3.3i_%2.2i.b',pthbin,expt,yr,iday,hr);


fld='thckn';
[F,n,m,l] = read_hycom(fina,finb,fld);
F(F>1e20)=nan;



fld='u-vel.';
[F,n,m,l] = read_hycom(fina,finb,fld);
F(F>1e20)=nan;
fn=1;
A=squeeze(F(1,:,:));
sub_plot_fld(HH,A,fn,fld);
plot(i0,j0,'k*');
contour(HH,[-3000:250:-100],'w')
caxis([-0.5 0.5]);

fld='v-vel.';
[F,n,m,l] = read_hycom(fina,finb,fld);
F(F>1e20)=nan;
fn=1;
A=squeeze(F(1,:,:));
sub_plot_fld(HH,A,fn,fld);
plot(i0,j0,'k*');
contour(HH,[-3000:250:-100],'w')
caxis([-0.5 0.5]);



