% Create frames
% Plot particles from Greenland Gyre
% see GreenlGyre_particles.m
% particles tracking 
%
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

slr = 1; % vertical layer
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
pthfig = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/%s/%3.3i/fig_GrGprtcl/',regn,expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/110/data_GG_prt/';
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

xlim1 = 200;
xlim2 = nn;
ylim1 = 5;
ylim2 = 1500;


YRPLT=[1993:1995];
nplt = length(YRPLT);


cc0=0;
% Combine all trajectories:
Xp=[];
Yp=[];
TM=[];
for ii=1:nplt
  yr = YRPLT(ii);
  fmat = sprintf('%sGG_particles-lr%2.2i_%i.mat',pthmat,slr,yr);
  fprintf('Loading saved %s\n',fmat);
  load(fmat);

  TR = PRTCL.TRACK;
  nr = length(TR);
  
  
  for it=1:nr    
    X = TR(it).I;
    Y = TR(it).J;
    Xp=[Xp,X];
    Yp=[Yp,Y];
    TM=[TM;TR(it).TM];
  end
  
end

TM=[datenum(1993,1,1);TM]; % missing initial date

[a1,a2]=size(Xp);
cmp=colormap_blue(a2);

Hmsk=HH*0;
Hmsk(HH<0)=1;

cc=0;
for itt=1:10:a2
  
  clf;
  pcolor(Hmsk); shading flat;
  colormap([0 0 0; 1 1 1]);
  hold on;
  
  contour(HH,[-8000:1000:-100],'Color',[0.5 0.5 0.5]);
  caxis([0 1]);
  
  plot(Xp(:,itt),Yp(:,itt),'.','Color',[0 0.5 1]);
  dv1=datevec(TM(itt));
%for ip=1:a1;
%  plot(Xp(ip,:),Yp(ip,:),'Color',cmp);
%  plot(Xp(ip,:),Yp(ip,:),'Color',[0 0.6 1]);
%end


  axis('equal');
  set(gca,'xlim',[xlim1 xlim2],...
	  'ylim',[ylim1 ylim2]);
  set(gca,'xtick',[],'ytick',[]);
    
  stt = sprintf('%2.2i/%2.2i/%i',dv1(3),dv1(2),dv1(1));
  title(stt);
    
  txtb = 'plot_GrGyre_prtcls.m';
  bottom_text(txtb,'pwd',1,'Fontsize',7);
%    drawnow
  cc = cc+1;
    
  if s_fig>0
    fgnm = sprintf('%s008-110_GrGyre_prtcl_%5.5i',pthfig,cc);
    fprintf('Saving figure %s\n',fgnm);
    print('-dpng','-r200',fgnm);
  end
  %fprintf('1 record: %8.5f min\n\n',toc/60);
		     %    keyboard

end

%keyboard  
  






