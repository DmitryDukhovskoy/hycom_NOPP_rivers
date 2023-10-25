% Plot trajectories of all partilces
% by layers
% all output have to be preprocessed in
% plot_AllGrSh_prtcls.m
%
% 
addpath /Net/Movies0/ddmitry/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /home/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers;
%addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
%addpath /Net/Movies0/ddmitry/MyMatlab/seawater
startup

close all
clear


SLR=[10,15,23,31]; % depth levels of the particels
ZLR=[50,90,150,450];  % nominal depths of particles in deep ocean
nlr=length(SLR);

%slr = 15; % vertical layer
%s_par=1; % parallel session 
%f_plt=1; % plot prtcles
s_fig = 0;
%cc = 71; % last SAVED frame - if need to restart from frame # 
      % the next frame will be cc+1
      % = 0 - start from beginning

% blkdat.input: target densities for the layers:
% upper 14 layers - fixed z-layers over the deep ocean
% 26.00   'sigma ' = layer 10 isopycnal target density (sigma units) - lr 5
% 30.65   'sigma ' = layer  A isopycnal target density (sigma units) - lr 15
% 35.20   'sigma ' = layer  I isopycnal target density (sigma units) - lr 23
% 36.70   'sigma ' = layer 22 isopycnal target density (sigma units) - lr 31

regn = 'ARCc0.08';
expt = 110;  
%pthfig  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/fig_trac/',expt);
pthfig = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/%s/%3.3i/fig_GrGprtcl/',regn,expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/110/data_GG_prt/';
txtb = 'plot_probability_prtcls.m';

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

xlim1 = 220;
xlim2 = 1250;
ylim1 = 30;
ylim2 = 1100;



fnmout=sprintf('%sGrShprt_comb.mat',pthmat);

fprintf('Loading %s\n',fnmout);
load(fnmout);

% Combine all experiments with Lagrangian particles together 
% for missing dates - use the previous number and locations of particles
nll=length(PP);
ntm=0;
for iip=1:nll
  tm=PP(iip).TM;
  nt=length(tm);
  ntm=max([ntm,nt]);
  if ntm==nt
    Li=iip;
  end
end

TM0=PP(Li).TM;
% Exclude repeated dates:
dTm=diff(TM0);
I=find(dTm>0);
I=[I;ntm];
TM0=TM0(I);
ntm=length(TM0);

% Fill in missing dates in the experiments 
% so that all experiments have same # of days (ntm)
clear nnm
nll=length(PP);
for iip=1:nll
  Xp=PP(iip).Xp;
  Yp=PP(iip).Yp;
  Tm=PP(iip).TM;
  zl=PP(iip).layer;
  [a1,a2]=size(Xp);

  Xn=[];
  Yn=[];
  for it=1:ntm
    t0=TM0(it);
    jt=find(Tm==t0,1);

    if ~isempty(jt),
      Xn(:,it)=Xp(:,jt);
      X0=Xp(:,jt);
      Yn(:,it)=Yp(:,jt);
      Y0=Yp(:,jt);
    else
      Xn(:,it)=X0;
      Yn(:,it)=Y0;
    end
  end

  PP(iip).Xp=Xn;
  PP(iip).Yp=Yn;
  PP(iip).ZL=ones(a1,1)*zl;
end


% Combine all experiments      
XP=[];
YP=[];
ZL=[]; 
TM=TM0;

for iip=1:nll
  Xp=PP(iip).Xp;
  Yp=PP(iip).Yp;
  zl=PP(iip).ZL;
  XP=[XP;Xp];
  YP=[YP;Yp];
  ZL=[ZL;zl];
end;




%[a1,a2]=size(XP);
%cmp=colormap_blue(a2);

% Colorcode for floats at different depths:
%CLRZ=[0 0.5 0.8;...
%     0.9 0.4 0;...
%     0. 0.8 0.2];  

Hmsk=HH*0;
Hmsk(HH<0)=1;

% Region of interest - double check
% with dS_FWC_timeseries_SubpolarGyre.m
IGR = [  430         729
         461         671
         559         660
         774         676
         781         608
         839         533
         855         530
         930         501
        1031         445
        1099         410
        1065         350
        1087         288
        1106         204
        1118         148
        1080         145
         444         145
         416         246
         367         492
         374         591
         379         729
        430       729];


% Count statistics
% of particles in the grid boxes
IL=find(HH>=0);
%SMM=zeros(nlr,mm,nn);
%for ilr=1:nlr
%  SMM(ilr,IL)=nan;
%end
fprintf('Calculating probability of Lagr partc. for layers ...\n');
for ilr=1:nlr 
  fprintf('Layer %i\n',ilr);

  SMM = zeros(mm,nn);
  SMM(IL)=nan;

		lr=SLR(ilr);
		IZ=find(ZL==lr);

  xmm=XP(IZ,:);
  ymm=YP(IZ,:);
  [a1,a2]=size(xmm);

  for iz=1:a1
    x0=round(xmm(iz,:));
    y0=round(ymm(iz,:));
    x0=x0(:);
    y0=y0(:);
    In=find(~isnan(x0));
    x0=x0(In);
    y0=y0(In);
    II=sub2ind(size(HH),y0,x0);
    SMM(II)=SMM(II)+1;
  end

  sgmx=3;
  sgmy=sgmx;
  npnts=3*sgmx;
  ilm1=xlim1-10;
  ilm2=xlim2+10;
  jlm1=ylim1-10;
  jlm2=ylim2+10;
  SMM = sub_gauss_filter(SMM,sgmx,sgmy,npnts,ilm1,ilm2,jlm1,jlm2);
% Scale to have overall integral = 1
  SMM = SMM/(nansum(nansum(SMM)));

  PR(ilr).SMM=SMM;
end
  
CMP = create_colormap_WBYR(400,0,1);
cmp = CMP.colormap;

%
for ilr=1:nlr
  figure(ilr); clf;
  pcolor(Hmsk); shading flat;
  colormap([0 0 0; 1 1 1]);
  hold on;
  
  caxis([0 1]);
  freezeColors;

  dmm = PR(ilr).SMM;
  dmm(HH>=0)=nan;
  I0=find(dmm==0);
  dmm(I0)=nan;
  dmm=log(dmm);
  dmm(I0)=-80;
  
  pcolor(dmm); shading flat;
		colormap(cmp);
  caxis([-18 -11]);

  plot(IGR(:,1),IGR(:,2),'-','Linewidth',1.5,...
       'Color',[0.6 0. 0.]);
  contour(HH,[-1000 -1000],'Color',[0.6 0.6 0.6]);
  contour(HH,[-500 -500],'Linewidth',1.6,'Color',[0.4 0.4 0.4]);
  
  axis('equal');
  set(gca,'xlim',[xlim1 xlim2],...
	  'ylim',[ylim1 ylim2]);
  set(gca,'xtick',[],'ytick',[]);
    
  clb = colorbar;
  set(clb,'Position',[0.91 0.11 0.017 0.81],...
          'TickLength',0.025,...
          'Fontsize',14);

  stt = sprintf('Prob Lagr prtcl in grid cell, log2, %i m',ZLR(ilr));
  title(stt);
    

  bottom_text(txtb,'pwd',1);
    
  if s_fig>0
%    fgnm = sprintf('%s008-110_GrSh_prtcl_lr%i',pthfig,ilr);
    fprintf('Saving figure %s\n',fgnm);
    print('-dpng','-r200',fgnm);
  end
  %fprintf('1 record: %8.5f min\n\n',toc/60);
		     %    keyboard

end





