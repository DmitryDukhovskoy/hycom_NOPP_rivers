% Create frames
% Plot particles from Greenland Gyre
% see GreenlGyre_particles.m
% particles tracking 
%  Also added particles 
% on NWest Labrador Shelf
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

%shelf='GrSh';
shelf='WLabr';

%NSM=[1,2,3,4,5,6];
NSM=1;
nnm=length(NSM);

slr = 31;
%slr = 15; % vertical layer
%slr = 23;
%slr = 31;
%f_plt=1; % plot prtcles
s_fig = 0;
%cc = 71; % last SAVED frame - if need to restart from frame # 
      % the next frame will be cc+1
      % = 0 - start from beginning

rg = 9806; 

dpth = -0;
switch(slr),
 case(10)
  dpth=-50;
 case(15)
  dpth=-90;
 case(23)
  dpth=-150;
 case(31);
  dpth=-450;
end

%DPTH=[-50; -90; -150; -450]; % approximate depth of the HYCOM layers in deep ocean

regn = 'ARCc0.08';
expt = 110;  
%pthfig  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/fig_trac/',expt);
pthfig = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/%s/%3.3i/fig_GrGprtcl/',regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/110/data_GG_prt/';
if strncmp(shelf,'WLabr',5)
  pthmat  = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/110/data_WLabr_prt/';
end
txtb = 'plot_GrSh_prtcls.m';

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
xlim2 = 1300;
ylim1 = 10;
ylim2 = 1300;


YRPLT=[1993:2019];
nplt = length(YRPLT);


cc0=0;
% Combine all trajectories:
for it=1:nnm
  PP(it).Xp=[];
  PP(it).Yp=[];
  PP(it).TM(1)=datenum(1993,1,1); % missing initial date
end

for ii=1:nplt
  yr = YRPLT(ii);
%  fmat = sprintf('%sGG_particles-lr%2.2i_%i.mat',pthmat,slr,yr);
% Combine all simulations
  for iip=1:nnm
    nsim=NSM(iip);
%    fmat = sprintf('%sGrSh_part-lr%2.2i_%i-%2.2i.mat',pthmat,slr,yr,nsim);
    fmat = sprintf('%s%s_part-lr%2.2i_%i-%2.2i.mat',pthmat,shelf,slr,yr,nsim);
    if exist(fmat,'file')
      fprintf('Loading saved %s\n',fmat);
      load(fmat);
    else
      fprintf('Does not exist %s\n',fmat);
    end
    
    
    TR = PRTCL.TRACK;
    nr = length(TR);
  
  
    Xp=PP(iip).Xp;
    Yp=PP(iip).Yp;
    TM=PP(iip).TM;
    for it=1:nr    
      X = TR(it).I;
      Y = TR(it).J;
      Xp=[Xp,X];
      Yp=[Yp,Y];
      TM=[TM;TR(it).TM];

    end
    PP(iip).Xp=Xp;
    PP(iip).Yp=Yp;
    PP(iip).TM=TM;
  end  % simulations
  
end


% Combine all experiments with Lagrangian particles together 
% for missing dates - use the previous number and locations of particles
ntm=0;
for iip=1:nnm
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
for iip=1:nnm
  Xp=PP(iip).Xp;
  Yp=PP(iip).Yp;
  Tm=PP(iip).TM;

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
end


% Combine all experiments      
XP=[];
YP=[];
TM=TM0;

for iip=1:nnm
  Xp=PP(iip).Xp;
  Yp=PP(iip).Yp;
  XP=[XP;Xp];
  YP=[YP;Yp];
end;




[a1,a2]=size(XP);
cmp=colormap_blue(a2);

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

%keyboard  
% Analyze # of particles in the SPG by years
% SPG region is defined in FWC_regions_GreenlExp.m
[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2

% N Atl regions
%BX = sub_define_boxes(HH,LON,LAT,1);
[XX,YY] = meshgrid((1:nn),(1:mm));


%INP = inpolygon(XX,YY,IGR(:,1),IGR(:,2));
%INGr= find(INP==1 & HH<0);
%INdeep = find(INP==1 & HH<-800);
%II = find(HH<0);

cc=0;
for itt=a2:a2
  figure(1); clf;
  pcolor(Hmsk); shading flat;
  colormap([0 0 0; 1 1 1]);
  hold on;
  
  contour(HH,[-8000:1000:-100],'Color',[0.8 0.8 0.8]);
  contour(HH,[-800 -800],'Linewidth',1.6,'Color',[0.5 0.5 0.5]);
  caxis([0 1]);
  
  plot(XP(:,itt),YP(:,itt),'.',...
       'Markersize',14,...
       'Color',[1 0.2 0]);
  
  plot(IGR(:,1),IGR(:,2),'-','Linewidth',2,...
       'Color',[0 0.8 0.4]);
  
  dv1=datevec(TM(itt));
%for ip=1:a1;
%  plot(Xp(ip,:),Yp(ip,:),'Color',cmp);
%  plot(Xp(ip,:),Yp(ip,:),'Color',[0 0.6 1]);
%end


  axis('equal');
  set(gca,'xlim',[xlim1 xlim2],...
	  'ylim',[ylim1 ylim2]);
  set(gca,'xtick',[],'ytick',[]);
    
  stt = sprintf('%2.2i/%2.2i/%i, Lr %i (%i m)',dv1(3),dv1(2),dv1(1),slr,dpth);
  title(stt);
    
  bottom_text(txtb,'pwd',1,'Fontsize',7);
%    drawnow
  cc = cc+1;
    
  if s_fig>0
    fgnm = sprintf('%s008-110_GrSh_prtcl_%5.5i',pthfig,cc);
    fprintf('Saving figure %s\n',fgnm);
    print('-dpng','-r200',fgnm);
  end
  %fprintf('1 record: %8.5f min\n\n',toc/60);
		     %    keyboard

end



TMd=TM-TM(1);
TMy=TMd./365.24;

clear INdp
Np0=length(find(~isnan(XP(:,1)))); % total # of particles initial
for ipp=1:a2
  NP(ipp,1)=length(find(~isnan(XP(:,ipp)))); % active particles total
% Particles in the domain:
  X0=XP(:,ipp);
  Y0=YP(:,ipp);
  II=inpolygon(X0,Y0,IGR(:,1),IGR(:,2));
  i0=length(find(II==1));
  INspg(ipp,1)=i0;
%
% Find particles in the deep basin
  ipt=round(X0(II==1));
  jpt=round(Y0(II==1));
  IJpt=sub2ind(size(HH),jpt,ipt);
  imm=find(HH(IJpt)<-800);
  INdp(ipp,1)=length(imm);
end

figure(2); clf;
axes('Position',[0.08 0.5 0.9 0.4]);
plot(TMy,INspg/Np0,'Linewidth',2);
hold
plot(TMy,INdp/Np0,'Linewidth',2);
stll=sprintf('Fraction of Lagr. prt in SPG and deep SPG, Lr %i (%i m)',slr,dpth);
title(stll);
set(gca,'Tickdir','out',...
	'Fontsize',14,...
	'xlim',[0 max(TMy)],...
	'xtick',[0:1:100],...
	'ylim',[0 1.05],...
	'ytick',[0:0.1:1.2],...
	'xgrid','on',...
	'ygrid','on');

xlabel('Years');
ylabel('Fraction');

bottom_text(txtb,'pwd',1,'Position',[0.08 0.3 0.4 0.04]);


% 
% Histogram of time inside the SPG by each particles
% Residnce time = average time a particles spends in the domain
npp=size(XP,1);
dTMy=diff(TMy);
dTMy(end+1)=dTMy(end);
for ipp=1:a1
  X0=XP(ipp,:);
  Y0=YP(ipp,:);
  II=inpolygon(X0,Y0,IGR(:,1),IGR(:,2));
  dtm=dTMy.*II';
  TMprt(ipp)=nansum(dtm);
end

xb=[1:2:26];
NN=hist(TMprt,xb);
NN=NN/a1; % fraction of all particles


figure(3); clf;
axes('Position',[0.08 0.5 0.9 0.4]);
hbr=bar(xb,NN,0.98);
set(hbr,'FaceColor',[0.6 0.6 0.6]);
stl2=sprintf('Time Lagr.prt in SPG, Lr %i (%im)',slr,dpth);
title(stl2);
set(gca,'Tickdir','out',...
	'Fontsize',14,...
	'xlim',[0 26],...
	'xtick',xb,...
	'ylim',[0 1.05*max(NN)],...
	'ytick',[0:0.05:1.2],...
	'xgrid','on',...
	'ygrid','on');

xlabel('Years');
ylabel('Fraction');

bottom_text(txtb,'pwd',1,'Position',[0.08 0.3 0.4 0.04]);





