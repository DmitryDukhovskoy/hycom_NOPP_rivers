% Read P. Myers's data
% of Fram fluxes
%
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/Wavelet;
startup;

close all
clear


YR1=2006;
YR2=2011;

Sref=34.8; 
rg  = 9806;
hgg = 1e20; % 

pthmat  = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/data_Myers/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';

yr=2018;
mo=12;
fnm=sprintf('%sANHA12-EXH006_y%4.4im%2.2id31_FRAMclip_gridT.nc',pthmat,yr,mo);

Z=-double(squeeze(nc_varget(fnm,'deptht')));
S=double(squeeze(nc_varget(fnm,'vosaline')));

ss=squeeze(S(:,1,:));
[a1,a2]=size(ss);
X=[1:a2];
for i=1:a2
  ib=max(find(~isnan(ss(:,i))));
  if ~isempty(ib)
    zb=Z(ib+1);
  else
    zb=0;
  end
  
  Hb(i)=zb;
end


figure(1); clf;
pcolor(X,Z,ss); shading flat;

stl=sprintf('arc08-%3.3i %s, S, %4.4i/%2.2i/%2.2i',...
	      expt,sctnm,dv(1:3));
stl= ' ';
cntr=[31:0.5:35];
cc0=34.8;
nfg=3;
mbtm=1;
sub_plot_Sxsct(nfg,Hb,ss,ZZ,X,dnmb,stl,cntr,cc0,mbtm);
