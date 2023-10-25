% Plot T or S fields
% Plot only 1 layer at a time
% if need to layer-average - modify the code
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

plr   = 1;  % layer to plot
pfld  = 'temp';
%pfld  = 'salin';
s_fig = 0;
rg = 9806;
fprintf('Plotting field: %s\n',pfld);

switch(pfld),
 case('salin');
%  c1=25.5;
%  c2=35.5;
   c1=29.5;
   c2=35.5;
%  c1=31;
%  c2=35.;
  ps=[1062 290 780 816];
 case('temp');
  c1=-2;
  c2=12;
  ps=[100 290 780 816];
end;


YRPLT=[];
cc=0;
for iyr=2017:2017
  for idd=245:245
    if idd==1, idd=2; end;
    cc=cc+1;
    YRPLT(cc,1)=iyr;
    YRPLT(cc,2)=idd;
  end
end

np=size(YRPLT,1);

% Experiments:
% 102 - test simulation, HYCOM GLBb0.08 nest: 1 file kept constant (1993,1,1)

regn = 'ARCc0.08';
%expt = 112;  % old GOFS3.1
expt = 123; % new sim with GOFS3.5

fprintf('%s_%3.3i\n',regn,expt);

%pthfig  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/fig_trac/',expt);
pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc0.08/%3.3i/data_mat/',expt);

inc1=1;
inc2=1600;
jnc1=1;
jnc2=2520;
djnc=2520;
dinc=1600;

ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2

xlim1 = 0;
xlim2 = nn;
ylim1 = 0;
ylim2 = mm;

% Plot fields:
cnc=0;
ip1=1;
for ip=ip1:np
  yr=YRPLT(ip,1);
  iday=YRPLT(ip,2);
  pthbin = sprintf('/nexsan/archive/%s_%3.3i/data/%i/',regn,expt,yr);
  if expt==112
    pthbin = sprintf('/nexsan/hycom/ARCc0.08_112/data/%i/',yr);
    fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
    finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);
  elseif expt==123
    pthbin='/nexsan/people/ddmitry/hycom/ARCc0.08_123/data/2017_dltEddTmltEAPJRA/';
    fina = sprintf('%s%3.3i_archv.%4.4i_%3.3i_00.a',pthbin,expt,yr,iday);
    finb = sprintf('%s%3.3i_archv.%4.4i_%3.3i_00.b',pthbin,expt,yr,iday);
  end

  
  if ~exist(fina,'file');
    fprintf('Not found: %s\n\n',fina);
    continue;
  end
  
  cnc=cnc+1;
  dnmb=datenum(yr,1,1)+iday-1;
  DV=datevec(dnmb);

  fprintf('%s, sfig=%i, V.Lr=%i, %4.4i_%2.2i_%2.2i: %s\n',...
	  pfld,s_fig,plr,DV(1:3),fina);
% in archm u_vel=utot; total velocity
%  keyboard
  f_getu=0;
  if f_getu==1
    ilr=1;
    [Fu,n,m,l] = read_hycom(fina,finb,'u-vel.');
    [Fv,n,m,l] = read_hycom(fina,finb,'v-vel.');
    U=squeeze(Fu(ilr,:,:));
    U(U>1e6)=nan;
    V=squeeze(Fv(ilr,:,:));
    V(V>1e6)=nan;
    S=sqrt(U.^2+V.^2);
  
  end

% Layer thickness:
  [F,n,m,l] = read_hycom(fina,finb,'thknss','r_layer',plr);
  F=squeeze(F);
  F=F./rg;
  F(F>1e10)=0;
  F(F<1e-2)=0;
  dH=squeeze(F); 

  tic;
  [F,n,m,l] = read_hycom(fina,finb,pfld,'r_layer',plr);
  toc;
  F(F>1e6)=nan;

%  F = squeeze(F);
%  lTr = log(F);
  lTr = squeeze(F);
  
  
  nf = 1;
  ifx=max(strfind(fina,'/'));
  zLdp = mean_ZM_41lrs(plr);
  lTr(HH>zLdp) = nan;
  stl=sprintf('%s, %s, Lr %i, %6.1fm %i/%2.2i/%2.2i',...
	      fina(ifx+1:end),pfld,plr,zLdp,DV(1:3));

% Greenland limits:
%  xlim1 = 400;
%  xlim2 = 1000;
%  ylim1 = 300;
%  ylim2 = 1100;  
% SW Greenland
  xlim1 = 375;
  xlim2 = 650;
  ylim1 = 300;
  ylim2 = 725;

  if strncmp(pfld,'salin',4),
    f_cmp=6;
  else
    f_cmp=3;
  end

  Fpos = [1730         485         754         849];

  close all
  hps = [0.9 0.1 0.02 0.8];
  sub_plot_scalar(lTr,nf,HH,xlim1,xlim2,...
      ylim1,ylim2,LON,LAT,stl,pfld,...
      'c1',c1,'c2',c2,'cmp',f_cmp,'clbpos',hps,'figpos',Fpos);


  txtb='plot_arc08_scalar.m';
  bottom_text(txtb,'pwd',1,'Position',[0.02 0.05 0.6 0.1]);
  
  set(gcf,'position',[1190 215 1138 1099]);
  
  if s_fig>0
    ffg=sprintf('%s%s',pthfig,fnmF);
    fprintf('Saving %s\n\n',ffg);
    print('-dpng','-r200',ffg);
%    sso=sprintf('./trim_jpeg.com %s %s',pthfig,fnmF);
%    system(sso);
  end
end;  % day loop

%I=find(~isnan(Tr));
%T=Tr(I);
%b=[1e-30;1e-10;1e-6;1e-3;1e-2;1e-1;1;(10:100)'];
%N=hist(T,b);
%bar(log10(b),N);




