% Plot difference of 
% Extract S fields and calc mean over specified days for given year 
% Monthly mean T/S NOT-averaged within the layers
% at the closest interface depths
% for 2 experiments (with & without Greenland runoff)
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

%lr1 = [0,-50];
%lr2 = [-50,-150];
%lr3 = [-150,-300]; % full water column

LRS = [0; -52; -152; -300; -500; -1000];
nlr = length(LRS)-1;
       

regn = 'ARCc0.08';
%expt = 110; % no Greenland runoff  
expt = 112;  % Greenland runoff



s_mat = 1; % =2 - load saved
pfld  = 'salin';

s_fig = 0;
rg = 9806;
fprintf('Extracting field: %s, Layer: %i\n',pfld);

% Averaging window
YRPLT=[];
cc=0;
for iyr=2014:2014  % do 1 year at a time
  for idd=330:5:365
    if idd==1, idd=2; end;
    cc=cc+1;
    YRPLT(cc,1)=iyr;
    YRPLT(cc,2)=idd;
  end
end
yr = YRPLT(1,1);
np=size(YRPLT,1);

% Experiments:
% 102 - test simulation, HYCOM GLBb0.08 nest: 1 file kept constant (1993,1,1)

pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
%pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc0.04/%3.3i/data_mat/',expt);

%fmat = sprintf('%sarc08_%3.3i_mnthly_%s_%i.mat',pthmat,expt,pfld,iyr);

%ftopo = sprintf('%s/depth_%s_17DD.nc',pthtopo,regn); % 
ftopo = sprintf('%s/depth_%s_09.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

%[DX,DY]=sub_dx_dy(LON,LAT);
%Acell=DX.*DY; % Grid cell area, m2


if s_mat==1
%  cnc=0;
  ip1=1;
  mold=0;

  for ilr=1:nlr
    S1(ilr).Fld=zeros(mm,nn);
    S1(ilr).nrec=0;
    S2(ilr).Fld=zeros(mm,nn);
    S2(ilr).nrec=0;
  end    

  for ip=ip1:np
    yr=YRPLT(ip,1);
    iday=YRPLT(ip,2);

    for ixp=1:2
      if ixp==1
	expt=110;
	TSFLD=S1;
	cnc=S1(1).nrec;
      else
	expt=112;
	TSFLD=S2;
	cnc=S2(1).nrec;
      end
      
      pthbin = sprintf('/nexsan/archive/%s_%3.3i/data/%i/',regn,expt,yr);
      if expt~=110,
	pthbin = sprintf('/nexsan/hycom/%s_%3.3i/data/%i/',regn,expt,yr);
      end

    %  pthbin = sprintf('/nexsan/hycom/ARCc0.04_011/data012/%i/',yr);  % Greenland on exp

      fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
      finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);

      if ~exist(fina,'file');
	fprintf('Not found: %s\n\n',fina);
	continue;
      end

      dnmb=datenum(yr,1,1)+iday-1;
      DV=datevec(dnmb);
      mo=DV(2);
      mday=DV(3);

      if mold==0; mold=mo; end;

      fprintf('%s, sfig=%i, %4.4i_%2.2i_%2.2i: %s\n',...
	      pfld,s_fig,DV(1:3),fina);

      tic;
%      cnc=cnc+1;
     % Layer thickness:
      [F,nn,mm,ll] = read_hycom(fina,finb,'thknss');
      F=squeeze(F);
      F=F./rg;
      F(F>1e10)=0;
      F(F<1e-2)=0;
      Lthck = F;

    % Create Depth array of interface depths:
    % Note these are BOTTOM interfaces 
    % So Layer 1 is between interfaces 0m and ZZ(1)
      ZZb = F.*nan;
      ZZb(1,:,:) = 0;
      I = find(HH>=0);
      ZZb(1,I) = nan;
      for kk=1:ll
	ZZb(kk+1,:,:)=ZZb(kk,:,:)-Lthck(kk,:,:);
      end

    % Read T/S:    
      [F,n,m,l] = read_hycom(fina,finb,pfld);
      F(F>1e6)=nan;
      Ctr = squeeze(F);
      clear F Lthck

    % Integrate tracer over the depths:
    % then divide by layer depths = mean tr. conc within 
    % specified depth levels
      fprintf('Integrating %s over depths ...\n',pfld);
      for ilr=1:nlr
	FLD(ilr,:,:)=zeros(mm,nn);
      end

      for k=1:ll
	aa0= squeeze(ZZb(k,:,:));
	aa = squeeze(ZZb(k+1,:,:));
	dZ = squeeze(abs(ZZb(k+1,:,:)-ZZb(k,:,:)));
    %    Inan = find(HH>=zbt); %lr(1)
	fm0 = squeeze(Ctr(k,:,:));
    %    fm0(Inan) = nan;

    % Go by layers:
	for ilr=1:nlr
	  zbt=LRS(ilr+1);
	  zup=LRS(ilr);
    % Interpolate onto depth level:
    % No interpolateion, S is within the layer
    % Find the layer that contains zbt
	  I = find(aa0>=zbt & aa<=zbt);
	  if isempty(I), continue; end;
  %	if ilr==4, keyboard; end;
    %      fmm = squeeze(Ctr(k+1,:,:));
    %      fmm(Inan) = nan;
    %      I = find(aa0>=zbt & aa<=zbt);
    %      dzL = zbt-aa0;
    %      fmI = (fmm-fm0)./(aa-aa0).*dzL+fm0;
	  fmI = squeeze(FLD(ilr,:,:));
	  fmI(I) = fm0(I);
	  FLD(ilr,:,:)=fmI;
%	  if ilr==4,
%	    fprintf('=== Lr=%i, k=%i, z=%5.1f, Check S= %8.6f\n',...
%		    ilr,k,aa(428,884),fmI(428,884));
%	  end
	  
	end  % LAYERS

      end      % vert layers

      for ilr=1:nlr
	cnc = TSFLD(ilr).nrec;
	fmI = squeeze(FLD(ilr,:,:));
	fprintf('=== Lr=%i    Max S= %8.6f\n',ilr,max(max(fmI)));
	fprintf('=== Lr=%i    Check S= %8.6f\n',ilr,fmI(428,884));

	TSFLD(ilr).TM       = dnmb;
	TSFLD(ilr).depth_av = [LRS(ilr), LRS(ilr+1)];
	TSFLD(ilr).Fld(:,:) = squeeze(TSFLD(ilr).Fld(:,:))+fmI;
	TSFLD(ilr).nrec     = cnc+1;
      end

      fprintf('Reading 1 rec: %8.6f min\n',toc/60);

      if ixp==1
	S1=TSFLD;
      else
	S2=TSFLD;
      end
    end; % expt loop
    
% Calculate differences:   
    f_chck=1;
    if f_chck==1
      fav=0;
      sub_anls_diff(S1,S2,HH,fav,mm,nn,nlr,LRS);
    end
  
  end;  % day loop

  for ilr=1:nlr
    cnc = TSFLD(ilr).nrec;
    dmm=TSFLD(ilr).Fld;
    TSFLD(ilr).Fld = dmm./cnc;
  end

  if s_mat==1
    fmat=sprintf('%sarc08_%3.3i_Sfld_zlv_%3.3i%3.3i.mat',pthmat,expt1,expt2,yr);
    fprintf('Saving %s\n',fmat);
    save(fmat,'S1','S2');
  end

end
  

%expt1=110;
%pthmat1  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt1);
%fmat1=sprintf('%sarc08_%3.3i_Sfld_zlv_%i.mat',pthmat1,expt1,yr);
%fprintf('Loading %s\n',fmat1);
%load(fmat1);
%S1=TSFLD;

%expt2=112;
%pthmat2  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt2);
%fmat2=sprintf('%sarc08_%3.3i_Sfld_zlv_%i.mat',pthmat2,expt2,yr);
%fprintf('Loading %s\n',fmat2);
%load(fmat2);
%S2=TSFLD;


cl1=colormap_blue(100);
cl2=colormap_red(100);
for ik=1:12;
  cl1(ik,:)=[1 1 1];
  cl2(ik,:)=[1 1 1];
end
cl1=flipud(cl1);
cmp = [cl1;cl2];
cmp = smooth_colormap(cmp,10);

c1 = -0.1;
c2 = 0.1;


xl1=300;
xl2=1200;
yl1=200;
yl2=1100;

Lmsk = HH*0;
Lmsk(HH<0)=1;
lmp=[0 0 0; 1 1 1];

for ilr=1:nlr
  F1 = S1(ilr).Fld;
  F2 = S2(ilr).Fld;
  F1(F1==0)=nan;
  F2(F2==0)=nan;
  dF = F2-F1;
  zb2 = LRS(ilr+1);
  dF(HH>zb2)=nan;
  
  figure(ilr); clf;
  pcolor(Lmsk); shading flat;
  colormap(lmp);
  freezeColors;
  hold on;
  
  pcolor(dF); shading flat;
  caxis([c1 c2]);
%  hold on;
%  contour(HH,[0 0],'w','linewidth',1.2);
  axis('equal');
  set(gca,'xlim',[xl1 xl2],...
	  'ylim',[yl1 yl2],...
	  'Color',[0.8 0.8 0.8],...
	  'xtick',[],...
	  'ytick',[]);
  colormap(cmp);
  chb=colorbar('location','eastoutside');
  set(chb,'Fontsize',14,'ticklength',0.02);
  
  ctl=sprintf('Fld: %s, dF=Gr-NoGr, lr=%i %i, %i',pfld,S1(ilr).depth_av(1),...
	      S1(ilr).depth_av(2),yr);
  title(ctl);
  
  btx='dltS_GreenldExp.m';
  bottom_text(btx,'pwd',1);
end


