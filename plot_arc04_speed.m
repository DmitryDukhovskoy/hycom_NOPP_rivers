% Plot speed  from ARCc0.04 
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

%plr   = 1;  % layer to plot
plr   = 16;  % layer to plot, ~100m - similar to Claudia
pfld  = 'speed';

s_fig = 0;
rg = 9806;
fprintf('Plotting field: %s\n',pfld);

YRPLT=[];
cc=0;
for iyr=2006:2006
  for idd=15:15
    if idd==1, idd=2; end;
    cc=cc+1;
    YRPLT(cc,1)=iyr;
    YRPLT(cc,2)=idd;
  end
end

np=size(YRPLT,1);

% Experiments:
% 102 - test simulation, HYCOM GLBb0.08 nest: 1 file kept constant (1993,1,1)

regn = 'ARCc0.04';
expt = 011;  
pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';
%pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc0.04/%3.3i/data_mat/',expt);

inc1=1;
inc2=1600;
jnc1=1;
jnc2=2520;
djnc=2520;
dinc=1600;

ftopo = sprintf('%s/depth_%s_17DD.nc',pthtopo,regn); % 
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
%  pthbin = '/Net/mars/ddmitry/hycom/ARCc0.04/output/';
  pthbin = sprintf('/nexsan/hycom/ARCc0.04_011/data/%i/',yr);  
  
  fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
  finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);
  
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
%    ilr=1;
  [Fu,n,m,l] = read_hycom(fina,finb,'u-vel.');
  [Fv,n,m,l] = read_hycom(fina,finb,'v-vel.');
  U=squeeze(Fu(plr,:,:));
  U(U>1e6)=nan;
  V=squeeze(Fv(plr,:,:));
  V(V>1e6)=nan;
  lTr=sqrt(U.^2+V.^2);


 % Layer thickness:
  [F,n,m,l] = read_hycom(fina,finb,'thknss','r_layer',plr);
  F=squeeze(F);
  F=F./rg;
  F(F>1e10)=0;
  F(F<1e-2)=0;
  dH=squeeze(F); 

  nf = 1;
  ifx=max(strfind(fina,'/'));
  zLdp = mean_ZM_41lrs(plr);
  lTr(HH>zLdp) = nan; 
  stl=sprintf('%s, %s, Lr %i, %6.1fm %i/%2.2i/%2.2i',...
	      fina(ifx+1:end),pfld,plr,zLdp,DV(1:3));

  sub_plot_scalar(lTr,nf,HH,xlim1,xlim2,ylim1,ylim2,LON,LAT,stl,pfld);

  txtb='plot_arc04_speed.m';
  bottom_text(txtb,'pwd',1,'Position',[0.2 0.05 0.7 0.1]);
  
  if s_fig>0
    fnmF=sprintf('trcr%2.2i_lr%2.2i_%i-%3.3i',nTr,plr,yr,iday);
    ffg=sprintf('%s%s',pthfig,fnmF);
    fprintf('Saving %s\n\n',ffg);
    print('-djpeg','-r200',ffg);
%    sso=sprintf('./trim_jpeg.com %s %s',pthfig,fnmF);
%    system(sso);
  end
end;  % day loop





