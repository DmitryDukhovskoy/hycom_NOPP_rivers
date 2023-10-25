% Combine all experiments with particles:
% N experiments with particles
% released at different depth levels
% 
% Compare different levels
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

fget=0;  % =0 - load previously loaded data

cc = 0;
ilr = 4;  % layer to plot

NSM=[1,2,3,4,5,6]; % each epxeriment released 100 particles at 1 depth
%NSM=[1,2,3,4];
nnm=length(NSM);

SLR=[10,15,23,31]; % depth levels of the particels
ZLR=[50,90,150,450];  % nominal depths of particles in deep ocean
nlr=length(SLR);

%YRPLT=[2013:5:2019];
YRPLT=1993;
nplt = length(YRPLT);


%slr = 15; % vertical layer
%s_par=1; % parallel session 
%f_plt=1; % plot prtcles

s_fig = 1;

fprintf('Save fig %i, Layer=%i\n',s_fig,ilr);

%cc = 71; % last SAVED frame - if need to restart from frame # 
      % the next frame will be cc+1
      % = 0 - start from beginning

rg = 9806; 

regn = 'ARCc0.08';
expt = 110;  
%pthfig  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/fig_trac/',expt);
pthfig = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/%s/%3.3i/fig_GrGprtcl/',regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/110/data_GG_prt/';
txtb = 'plot_GrShprt_paths.m';

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

% layer interface depths:
% in N Atl:
zz=[                     0
                        -1
         -2.80000007967061
         -6.04000002390118
         -10.7199998326917
         -15.6499996414823
         -21.4599995777458
         -28.3299994502728
         -36.3299994502728
         -44.3299994502728
         -52.3299994502728
         -60.3299994502728
         -68.3299994502728
         -76.3299994502728
         -84.3299994502728
         -92.3299994502728
         -100.329999450273
         -108.329999450273
         -116.329999450273
         -124.329999450273
         -132.329999450273
         -140.329999450273
         -148.329999450273
         -156.329999450273
         -166.329999450273
         -182.730000087638
         -218.650001234894
         -261.030001362367
         -311.050001872259
         -370.070002382151
         -439.709999577746
         -521.889997793124
         -642.060033995449
         -833.055521452108
         -1631.30343089531
         -2130.04884186818
          -2911.1784563899
         -3125.48785879659
         -3385.88090896996
         -3385.88090896996
         -3385.88090896996
         -3385.88090896996];

f_dps = 0;
if f_dps==1
% Read layer thicknesses:
  pthbin = '/nexsan/hycom/ARCc0.08_112/data/1993/';
  fina = sprintf('%s112_archm.1993_002_12.a',pthbin);
  finb = sprintf('%s112_archm.1993_002_12.b',pthbin);

  [ZM,ZZ] = sub_zz_zm(fina, finb,HH);
% on the shelf depth of the layer:
 ZM(31,539,529)
 
end


fnmout=sprintf('%sGrShprt_comb.mat',pthmat);

if fget>0
  % Combine all trajectories:
  ik=0;
  for ip=1:nlr
    for it=1:nnm
      ik=ik+1;
      PP(ik).Xp=[];
      PP(ik).Yp=[];
      PP(ik).TM(1)=datenum(1993,1,1); % missing initial date

    end
  end


  for ii=1:nplt
    yr = YRPLT(ii);
  %  fmat = sprintf('%sGG_particles-lr%2.2i_%i.mat',pthmat,slr,yr);
  % Combine all simulations
    icc=0;  % all nlevels*nexpts
    for ilr=1:nlr
      slr=SLR(ilr);
      for iip=1:nnm
	       nsim=NSM(iip);
       	fmat = sprintf('%sGrSh_part-lr%2.2i_%i-%2.2i.mat',pthmat,slr,yr,nsim);
	
								if ~exist(fmat,'file')
										fprintf('Not found %s\n',fmat);
										continue;
								end
								
								fprintf('Loading saved %s\n',fmat);
								load(fmat);


								TR = PRTCL.TRACK;
								nr = length(TR);

								icc=icc+1;
								
								Xp=PP(icc).Xp;
								Yp=PP(icc).Yp;
								TM=PP(icc).TM;
								for it=1:nr    
										X = TR(it).I;
										Y = TR(it).J;
										Xp=[Xp,X];
										Yp=[Yp,Y];
										TM=[TM;TR(it).TM];
								end
								PP(icc).Xp=Xp;
								PP(icc).Yp=Yp;
								PP(icc).TM=TM;
								PP(icc).layer=slr;
								PP(icc).nsim=nsim;
      end  % simulations in 1 layer
    end  % layers
  end

  fprintf('Saving %s\n',fnmout);
  save(fnmout,'PP');
  
else
  fprintf('Loading %s\n',fnmout);
  load(fnmout);
end


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




[a1,a2]=size(XP);
nn1=a1/4;         % # of partcles in 1 layer
cc1=1;
cc2=nn1;
CMP = create_colormap_paleB(nn1,cc1,cc2);
CLRF = CMP.colormap;

% Colorcode for floats at different depths:
%CLRZ=[0 0.5 0.8;...
%     0.9 0.4 0;...
%     0. 0.8 0.2;...
%     0.5 0.1 1];  

Hmsk=HH*0;
Hmsk(HH<0)=1;

% Region of interest - double check
% with dS_FWC_timeseries_SubpolarGyre.m
fspg='SPG_noNorth_indx.mat';  % get IGR indices of the region
load(fspg);
IGR(end+1,:)=IGR(1,:);


% Particles release location:
%fmatLP = sprintf('%sGrSh_part-lr10_1993-01.mat',pthmatLP);
%Ip = PRTCL.TRACK(1).I;
%Jp = PRTCL.TRACK(1).J;
slr=[];
fina=[];
finb=[];
Np0=1e4;
nsim=1;
[bmm,dmm,amm]=sub_seed_GrSh_prt(slr,HH,LON,LAT,fina,finb,Np0,nsim);
Ip = bmm.TRACK.I;
Jp = bmm.TRACK.J;

%
% Labrador Sea:
[bmm,dmm,amm]=sub_seed_WLabr_prt(slr,HH,LON,LAT,fina,finb,Np0,nsim);
Ilb = bmm.TRACK.I;
Jlb = bmm.TRACK.J;


%keyboard  
% Analyze # of particles in the SPNA by years
% SPNA region is defined in FWC_regions_GreenlExp.m
%[DX,DY]=sub_dx_dy(LON,LAT);
%Acell=DX.*DY; % Grid cell area, m2

% N Atl regions
%BX = sub_define_boxes(HH,LON,LAT,1);
%[XX,YY] = meshgrid((1:nn),(1:mm));


%INP = inpolygon(XX,YY,IGR(:,1),IGR(:,2));
%INGr= find(INP==1 & HH<0);
%INdeep = find(INP==1 & HH<-800);
%II = find(HH<0);

%cc=0;
figure(1);
set(gcf,'Position',[1203 388 905 954]);
%set(gcf,'Position',[9 5 707 797]);

DV   = datevec(TM);
nrc  = length(TM);
%for itt=733:733
for iyr=1:nplt
  yr0=YRPLT(iyr);
  I=find(DV(:,1)== yr0);
  it1=I(1);
  it2=I(end);
  itt=it2;
  dv1=datevec(TM(it2));
  fprintf('Plotting %i/%i/%i\n',dv1(1:3));
  
  figure(1); clf;
  pcolor(Hmsk); shading flat;
  colormap([0 0 0; 1 1 1]);
  freezeColors;
  caxis([0 1]);

  hold on;
  
  contour(HH,[-8000:1000:-100],'Linewidth',1,'Color',[0.8 0.8 0.8]);
  contour(HH,[-1000 -1000],'Linewidth',1.6,'Color',[0.5 0.5 0.5]);
  contour(HH,[-500 -500],'Linewidth',1,'Color',[0.8 0.8 0.8]);

% Plot Lagr Prt WLabr Shelf Particles
% release locations
  plot(Ip,Jp,'.','Color',[0.7 0.8 1]);
  plot(Ilb,Jlb,'.','Color',[1 0.8 0.7]);

 
%keyboard
		lr=SLR(ilr);
		IZ=find(ZL==lr);
% do not plot pathways
%		for izz=1:length(IZ);
%				if mod(izz,10)==0,
%						fprintf('Partcles %4.1f%% ...\n',izz/length(IZ)*100);
%				end
%				iz0=IZ(izz);
%%				it1=itt-183+1;
%%				it2=itt;
%				for itm=it1:it2
%						dv = DV(itm,:);
%						dj1 = datenum(dv(1),1,1);
%						jday = TM(itm)-dj1+1;
%						clr = CLRF(izz,:);
%						plot(XP(iz0,itm:itm+1),YP(iz0,itm:itm+1),'-','Color',clr);
%				end
%		end
%
% Plot particles
		plot(XP(IZ,itt),YP(IZ,itt),'.',...
					'Markersize',28,...
					'Color',[0 0.4 0.8]);

  plot(IGR(:,1),IGR(:,2),'-','Linewidth',1.6,...
       'Color',[0.6 0. 0.]);
  
%for ip=1:a1;
%  plot(Xp(ip,:),Yp(ip,:),'Color',cmp);
%  plot(Xp(ip,:),Yp(ip,:),'Color',[0 0.6 1]);
%end

  axis('equal');
  set(gca,'xlim',[xlim1 xlim2],...
	  'ylim',[ylim1 ylim2]);
  set(gca,'xtick',[],'ytick',[]);
    
  stt = sprintf('%2.2i/%2.2i/%i, Lr=%i, Z=%i m',dv1(3),dv1(2),dv1(1),SLR(ilr),ZLR(ilr));
  title(stt);
  
  
  bottom_text(txtb,'pwd',1,'Fontsize',7);
%    drawnow
  cc = cc+1;
    
  if s_fig>0
    fgnm = sprintf('%s008-110_GrShLr%i_pths_%4.4i',pthfig,ilr,cc);
    fprintf('Saving figure %s\n',fgnm);
    print('-dpng','-r200',fgnm);
  end
  %fprintf('1 record: %8.5f min\n\n',toc/60);
		     %    keyboard

end







