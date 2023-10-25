% Create frames
% Plot particles from BG
% see BG_particles.m
% particles tracking 
% seeded in the upper BG
% Particles are not added during the simulation
% All N prticles seeded at once at initial state
% For each particle: 
% no T or S is tracked, only time and location
% 
addpath /Net/Movies0/ddmitry/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /home/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers;
%addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
%addpath /Net/Movies0/ddmitry/MyMatlab/seawater
startup

close all
clear

%s_par=1; % parallel session 
%f_plt=1; % plot prtcles
s_fig = 1;
cc = 71; % last SAVED frame - if need to restart from frame # 
      % the next frame will be cc+1
      % = 0 - start from beginning

rg = 9806; 

regn = 'ARCc0.08';
expt = 110;  
%pthfig  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/fig_trac/',expt);
pthfig = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/%s/%3.3i/fig_BGprtcl/',regn,expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/110/data_BG_prt/';
% ------------------------
% TOPO
% ------------------------
ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

hmsk=HH;
hmsk(HH<0)=nan;

xlim1 = 20;
xlim2 = nn;
ylim1 = 5;
ylim2 = 2000;


YRPLT=[1994:2003];
nplt = length(YRPLT);

c1=30;
c2=35;

%close all
%figure(1); clf;
%ff=figure('visible','off');

cc0=0;
for ii=1:nplt
  yr = YRPLT(ii);
  fmat = sprintf('%sBG_particles_%i.mat',pthmat,yr);
  fprintf('\n   Loading saved %s\n\n',fmat);
  load(fmat);

  TR = PRTCL.TRACK;
  nr = length(TR);
  
  
  for it=1:nr
    cc0=cc0+1;
    if cc0<cc, 
      fprintf('Searching restart position %i, %3.3i\n',cc,cc0);
      continue; 
    end;
    
tic;
    dnmb = TR(it).TM;
    if isempty(dnmb),
      dnmb = datenum(yr,1,1);
    end
    
    dv1 = datevec(dnmb);
    yr = dv1(1);
    mo = dv1(2);
    md = dv1(3);
    iday = dnmb-datenum(yr,1,1)+1;
% Plot SSS     
    pthbin = sprintf('/nexsan/archive/ARCc0.08_110/data/%i/',yr)
    fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
    finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);
    
    [F,nn,mm,ll] = read_hycom(fina,finb,'salin','r_layer',1);
    F = squeeze(F);
    F(F>1e6)=nan;

    if ii==-1e9
      Xm1=TR(it).I;
      Ym1=TR(it).J;
      Xm2=Xm1;
      Ym2=Ym1;
      Xm3=Xm1;
      Ym3=Ym1;
    end
  
   close all
%figure(1); clf;
   ff=figure('visible','off');
   clf;
%    contour(HH,[0 0],'k');
    pcolor(F); shading flat;
    caxis([c1 c2]);
    hold on;
    

    X = TR(it).I;
    Y = TR(it).J;
    plot(X,Y,'m.');
    in=length(X);
%    Xm1 = Xm1(1:in);
%    Ym1 = Ym1(1:in);
%    Xm2 = Xm2(1:in);
%    Ym2 = Ym2(1:in);
%    Xm3 = Xm3(1:in);
%    Ym3 = Ym3(1:in);
%    xa=[X Xm1]';
%    ya=[Y Ym1]';
%    plot(xa,ya,'m-');
%    xa=[Xm1 Xm2]';
%    ya=[Ym1 Ym2]';
%    plot(xa,ya);
%    xa=[Xm2 Xm3]';
%    ya=[Ym2 Ym3]';
%    plot(xa,ya);
    
%    Xm3=Xm2;
%    Xm2=Xm1;
%    Xm1=X;
%    Ym3=Ym2;
%    Ym2=Ym1;
%    Ym1=Y;
%keyboard
    axis('equal');
    set(gca,'xlim',[xlim1 xlim2],...
	    'ylim',[ylim1 ylim2],...
	    'Color',[0. 0. 0.]);
    set(gca,'xtick',[],'ytick',[]);
    
    hb = colorbar;
    set(hb,'Position',[0.75 0.11 0.03 0.8],...
	   'TickLength',0.042);
    
    stt = sprintf('%2.2i/%2.2i/%i',dv1(3),dv1(2),dv1(1));
    title(stt);
    txtb = 'plot_BGprtcls.m';
    bottom_text(txtb,'pwd',1,'Fontsize',7);
%    drawnow
    cc = cc+1;
    
    if s_fig>0
      fgnm = sprintf('%sARCc008_BGprtcl_%5.5i',pthfig,cc);
      fprintf('Saving figure %s\n',fgnm);
      print('-dpng','-r250',fgnm);
    end
    fprintf('1 record: %8.5f min\n\n',toc/60);
%    keyboard
  end
%keyboard  
  
end





