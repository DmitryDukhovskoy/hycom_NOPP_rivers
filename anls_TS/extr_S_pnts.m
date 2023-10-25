% Extract S at some locations
% 
% Monthly mean T/S averaged within the layers
% for 2 experiments (with & without Greenland runoff)
% expt_110 - no Greenland runoff
% expt 112 - with Greenland runoff
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

fprintf('expt %3.3i\n\n',expt);


s_mat = 2; % =2 - load and start from last saved
pfld  = 'salin';

s_fig = 0;
rg = 9806;
fprintf('Extracting field: %s, Layer: %i\n',pfld);

% Averaging window
%YRPLT=[];
%cc=0;
%for iyr=1993:2016  % do 1 year at a time
%  for idd=1:5:365
%    if idd==1, idd=2; end;
%    cc=cc+1;
%    YRPLT(cc,1)=iyr;
%    YRPLT(cc,2)=idd;
%  end
%end
%np=size(YRPLT,1);

% Experiments:
% 102 - test simulation, HYCOM GLBb0.08 nest: 1 file kept constant (1993,1,1)

pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);

fmat = sprintf('%sarc08_%3.3i_%s_sections.mat',pthmat,expt,pfld);

ftopo = sprintf('%s/depth_%s_09.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

%[DX,DY]=sub_dx_dy(LON,LAT);
%Acell=DX.*DY; % Grid cell area, m2
clear SCT
SCT = sub_gsa_sections;
nst = length(SCT);

for ik=1:nst
  xy=SCT(ik).XY;
  IJ=[];
  for k=1:2
    x0=xy(k,1);
    y0=xy(k,2);
    D=distance_spheric_coord(LAT,LON,y0,x0);
    [j,i]=find(D==min(min(D)),1);
    IJ(k,1)=i;
    IJ(k,2)=j;
  end
  SCT(ik).IJ=IJ;
end

for ik=1:nst
  IJ=SCT(ik).IJ;
  i1=IJ(1,1);
  i2=IJ(2,1);
  j1=IJ(1,2);
  j2=IJ(2,2);
  [I,J]=sub_xsct_indx(i1,j1,i2,j2);
  SCT(ik).IIs=I;
  SCT(ik).JJs=J;
  SCT(ik).IND=sub2ind(size(HH),J,I);
end

% Plot
f_chck=0;
if f_chck==1
  figure(10); clf;
  contour(HH,[0 0],'k');
  hold on;
  contour(HH,[-500 -500],'b');
  contour(HH,[-10000:1000:-900],'c');
  axis('equal');
  
  for ik=1:nst
    II=SCT(ik).IIs;
    JJ=SCT(ik).JJs;
    plot(II,JJ,'r.-');
    text(II(1),JJ(1),sprintf('%2.2i',ik),'Fontsize',12);
  end
end

dv0=[1993,0,1];
cnc = 0;
if s_mat==2
  fprintf('Loading existing %s\n',fmat);
  load(fmat);
  TM=SCT(end).TM;
  cnc=length(TM);
  fprintf('Found last record #i: %s\n',cnc,datestr(TM(end)));
  dv0=datevec(TM(end));
  ddn=32-dv0(3);
  dv2=datevec(TM(end)+ddn);
  fprintf('Will continue extraction from %i/%2.2i\n\n',dv2(1:2));
end
YR1=dv0(1);
MM1=dv0(2);
%  keyboard
for iyr=YR1:2016
  yr=iyr;
  dd=7;

  pthbin = sprintf('/nexsan/archive/%s_%3.3i/data/%i/',regn,expt,yr);
  if expt~=110,
    pthbin = sprintf('/nexsan/hycom/%s_%3.3i/data/%i/',regn,expt,yr);
  end

  for im=1:12
%  yr=YRPLT(ip,1);
%  iday=YRPLT(ip,2);
    if (iyr==YR1 & im<=MM1), continue; end;
    dd1=datenum(yr,im,1);
    dd2=dd1+31;
    dv=datevec(dd2);
    dd2=datenum(dv(1),dv(2),1)-1;
    ndays = dd2-dd1+1;

    for ilr=1:nlr
      FLD(ilr,:,:)=zeros(mm,nn);
    end

    cc=0;
    for mday=1:dd:ndays
      dnmb = datenum(yr,im,mday);
      dv = datevec(dnmb);
      iday =dnmb-datenum(dv(1),1,1)+1;

      fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
      finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);

      if ~exist(fina,'file');
	fprintf('Not found: %s\n\n',fina);
	continue;
      end
      if ~exist(finb,'file');
	fprintf('Not found: %s\n\n',finb);
	continue;
      end

      dnmb=datenum(yr,1,1)+iday-1;
      DV=datevec(dnmb);
      mo=DV(2);
      mday=DV(3);

      fprintf('%s, sfig=%i, %4.4i_%2.2i_%2.2i: %s\n',...
	      pfld,s_fig,DV(1:3),fina);

      tic;
      cc=cc+1; % day counter
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

      for k=1:ll
	aa0= squeeze(ZZb(k,:,:));
	aa = squeeze(ZZb(k+1,:,:));
	dZ = squeeze(abs(ZZb(k+1,:,:)-ZZb(k,:,:)));
	fm0 = squeeze(Ctr(k,:,:));

    % Go by layers:
	for ilr=1:nlr
	  zbt=LRS(ilr+1);
	  zup=LRS(ilr);
	  I = find(aa0>=zbt & aa<=zbt);
	  if isempty(I), continue; end;
	  fmI = squeeze(FLD(ilr,:,:));
	  fmI(I) = fmI(I)+fm0(I);
	  FLD(ilr,:,:)=fmI;

	end  %specified   LAYERS
      end      % model vert layers

      fprintf('Reading 1 rec: %8.6f min\n\n',toc/60);
      
    end  % days in month
  
    cnc=cnc+1; % monthly means
    for ik=1:nst
      SCT(ik).TM(cnc)=dnmb;
    end
%
% Save sections
    for ik=1:nst
      SCT(ik).TM(cnc)=dnmb;
      II=SCT(ik).IND;
      for ilr=1:nlr
	fmI = squeeze(FLD(ilr,:,:))/cc;
	fmI(fmI==0)=nan;
	ss=nanmean(fmI(II));
%	fprintf('section %i, Layer=%i, S=%5.3f\n',ik, ilr,ss);
	SCT(ik).SS(ilr,cnc) = ss;
	smin = min(fmI(II));
	SCT(ik).Smin(ilr,cnc) = smin;
	fprintf('sect %i, Lr=%i, Sav=%5.3f, Smin=%5.3f\n',ik, ilr,ss,smin);
      end
    end


  end; % months

  if s_mat>0
    fprintf('Saving %s\n',fmat);
    save(fmat,'SCT');
  end
  
end;  % day loop


if s_mat>0
%  fmat=sprintf('%sarc08_%3.3i_%s_sections.mat',pthmat,expt,pfld);
  fprintf('Saving %s\n',fmat);
  save(fmat,'SCT');
end


  



