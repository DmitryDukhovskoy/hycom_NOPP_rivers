% HYCOM-CICEv5 GOFS3.5 system
% new simulations 
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

f_tiles = 1;

plr   = 1;  % layer to plot
%plr   = 16;  % layer to plot, ~100m - similar to Claudia
pfld  = 'speed';

s_fig = 0;
rg = 9806;
fprintf('Plotting field: %s\n',pfld);

YRPLT=[];
cc=0;
for iyr=2017:2017
  for idd=161:161
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
expt = 022;  
%pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.04/topo_grid/';
%pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc0.04/%3.3i/data_mat/',expt);

%inc1=1;
%inc2=1600;
%jnc1=1;
%jnc2=2520;
%djnc=2520;
%dinc=1600;

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


if f_tiles == 1
  [XT,YT] = read_tiles;
end

% Plot fields:
cnc=0;
ip1=1;
for ip=ip1:np
  yr=YRPLT(ip,1);
  iday=YRPLT(ip,2);
%  pthbin = sprintf('/nexsan/hycom/ARCc0.04_011/data/%i/',yr); 
%  pthbin = '/Net/kronos/ddmitry/hycom/ARCc0.04/output/'; 
  pthbin = sprintf('/nexsan/people/ddmitry/hycom/ARCc0.04_022/data/%i/',yr);
%   pthbin = '/Net/kronos/ddmitry/hycom/ARCc0.04/output_eloan/';  % Energy Loan
 
%  fina = sprintf('%sarchv.%4.4i_%3.3i_00.a',pthbin,yr,iday);
%  finb = sprintf('%sarchv.%4.4i_%3.3i_00.b',pthbin,yr,iday);
%  fina = sprintf('%sarchm.%4.4i_%3.3i_12.a',pthbin,yr,iday); % mean fields - saved U=u'+u_btrop, total
%  finb = sprintf('%sarchm.%4.4i_%3.3i_12.b',pthbin,yr,iday);

  fina = sprintf('%s022_archv.%4.4i_%3.3i_00.a',pthbin,yr,iday);
  finb = sprintf('%s022_archv.%4.4i_%3.3i_00.b',pthbin,yr,iday);

  
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
% in archv - need to add barotropic u
%  keyboard
%    ilr=1;
  [A,n,m,l] = read_hycom(fina,finb,'u_btrop');
  [B,n,m,l] = read_hycom(fina,finb,'v_btrop');
  Ub=squeeze(A);
  Ub(Ub>1e6)=nan;
  Vb=squeeze(B);
  Vb(Vb>1e6)=nan;

  [Fu,n,m,l] = read_hycom(fina,finb,'u-vel.','r_layer',plr);
  [Fv,n,m,l] = read_hycom(fina,finb,'v-vel.','r_layer',plr);
  U=squeeze(Fu);
  U(U>1e6)=nan;
  Ut=U+Ub;
  V=squeeze(Fv);
  V(V>1e6)=nan;
  Vt=V+Vb;

  SS=sqrt(Ut.^2+Vt.^2);


 % Layer thickness:
%  [F,n,m,l] = read_hycom(fina,finb,'thknss','r_layer',plr);
%  F=squeeze(F);
%  F=F./rg;
%  F(F>1e10)=0;
%  F(F<1e-2)=0;
%  dH=squeeze(F); 

  nf = 1;
  ifx=max(strfind(fina,'/'));
%  zLdp = mean_ZM_41lrs(plr);
%  SS(HH>zLdp) = nan; 
  stl=sprintf('%s, %s, Lr %i, %i/%2.2i/%2.2i',...
	      fina(ifx+1:end),pfld,plr,DV(1:3));

%  sub_plot_scalar(SS,nf,HH,xlim1,xlim2,ylim1,ylim2,LON,LAT,stl,pfld);

  txtb='plot_arc04_speed.m';
%  bottom_text(txtb,'pwd',1);

  f_bcl = 1;  % plot bacl U,V
  if f_bcl==1
    figure(5); clf;
    pcolor(U); shading flat;
    hold on;

    colorbar
    caxis([-0.5 0.5]);

    if f_tiles == 1 % plot tiles:
      NN=length(XT);
      ctt = [0.7 0.7 0.7];
      for ii=1:NN
        i1=XT(ii,1);
        i2=XT(ii,2);
        j1=YT(ii,1);
        j2=YT(ii,2);    
        plot([i1 i2],[j1 j1],'color',ctt);
        plot([i1 i2],[j2 j2],'color',ctt);
        plot([i1 i1],[j1 j2],'color',ctt);
        plot([i2 i2],[j1 j2],'color',ctt);
      end
    end

    stl=sprintf('%s, u-anom, %i %2.2i/%2.2i/%2.2i',...
        fina(ifx+1:end),plr,DV(1:3));
    title(stl,'Interpreter','None');
    bottom_text(txtb,'pwd',1);


    figure(6); clf;
    pcolor(V); shading flat;
    hold on
    colorbar
    caxis([-0.5 0.5]);

    if f_tiles == 1 % plot tiles:
      NN=length(XT);
      ctt = [0.7 0.7 0.7];
      for ii=1:NN
        i1=XT(ii,1);
        i2=XT(ii,2);
        j1=YT(ii,1);
        j2=YT(ii,2);    
        plot([i1 i2],[j1 j1],'color',ctt);
        plot([i1 i2],[j2 j2],'color',ctt);
        plot([i1 i1],[j1 j2],'color',ctt);
        plot([i2 i2],[j1 j2],'color',ctt);
      end
    end

    stl=sprintf('%s, v-anom, Lr %i %2.2i/%2.2i/%2.2i',...
        fina(ifx+1:end),plr,DV(1:3));
    title(stl,'Interpreter','None');
    bottom_text(txtb,'pwd',1);

  end

  f_btrop=0;
  if f_btrop==1

    figure(10); clf;
    pcolor(Ub); shading flat;
    colorbar
    caxis([-0.5 0.5]);

    stl=sprintf('%s, u-btrop, %2.2i/%2.2i',...
        fina(ifx+1:end),DV(1:3));

    bottom_text(txtb,'pwd',1);


    figure(11); clf;
    pcolor(Vb); shading flat;
    colorbar
    caxis([-0.5 0.5]);

    stl=sprintf('%s, v-btrop, %2.2i/%2.2i',...
        fina(ifx+1:end),DV(1:3));
    bottom_text(txtb,'pwd',1);


  end

  
end;  % day loop





