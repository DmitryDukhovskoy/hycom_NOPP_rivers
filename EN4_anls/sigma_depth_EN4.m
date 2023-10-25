% Plot deoths of sigma surfaces
% http://hadobs.metoffice.com/en4/index.html
%The EN4 dataset consists of two products:
%Observed subsurface ocean temperature and 
% salinity profiles with data quality information, and,
%Objective analyses formed from the profile data with uncertainty estimates.
%Data are available from 1900 to the present and 
% there are separate files for each month.
% Please read 'Good, S. A., M. J. Martin and N. A. Rayner, 2013. EN4: 
% quality controlled ocean temperature and salinity profiles and 
% monthly objective analyses with uncertainty estimates, 
% Journal of Geophysical Research: Oceans, 118, 6704-6716, 
% doi:10.1002/2013JC009067' for details of how the dataset was constructed.
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

s_fig = 1; 
sgm0  = 28;  % sigma surface 
s_mat = 1;
s_calc = 0; % derive depths

%pthdat = '/Net/yucatan/tachanat/ocean_analysis/EN4/EN4_extract/';
pthdat = '/Net/kronos/ddmitry/EN4/';
pthmat = '/Net/ocean/ddmitry/vector_winds/dataM/'; 
pthfig = '/Net/ocean/ddmitry/vector_winds/fig_EN4/';
en4v   = 'EN.4.2.0.f.analysis.g10'; % EN4 version
txtb   = 'anls_mtlb_utils/vector_winds/sigma_depth_EN4.m';
fmat   = sprintf('%sdepth_sigma_%4.4i_EN4.mat',pthmat,round(sgm0*100));

%YC1=1993;
%YC2=1993;
yr1=1993;
yr2=2015;
%fclim = sprintf('%sEN4_climatology%i_%i.mat',pthmat,YC1,YC2);
%fanom = sprintf('%sEN4_climatology%i_%i.mat',pthmat,yr1,yr2);

fnm = sprintf('%s%s.201302.nc',pthdat,en4v);

S = double(nc_varget(fnm,'salinity'));
S = squeeze(S(1,1,:,:));

% Find sections:
LON = nc_varget(fnm,'lon');
LAT = nc_varget(fnm,'lat');
ZZ  = -1*nc_varget(fnm,'depth');
I = find(LON>180);
LON(I)=LON(I)-360;

% Depth in the central Greenland Sea: 
i1=max(find(LON<=-10));
j1=min(find(LAT>=70));
i2=min(find(LON>=10));
j2=max(find(LAT<=76));

SGMD.Name   = sprintf('GrSea, Depth Sigma surf=%6.2f',sgm0);
SGMD.Code   = 'anls_mtlb_utils/vector_winds/sigma_depth_EN4.m';
SGMD.IJ     = [i1,j1;i2,j2];
SGMD.YCoord = [LAT(j1),LAT(j2)];
SGMD.XCoord = [LON(i1),LON(i2)];


cc=0;
for ii=yr1:yr2
  for im=1:12
    cc=cc+1;
    YRPLT(cc,1)=ii;
    YRPLT(cc,2)=im;
  end
end

if s_calc>0
  nrc=cc;
  cc=0;
  for ik=1:nrc
    YR=YRPLT(ik,1);
    IM=YRPLT(ik,2);
    fnm = sprintf('%s%s.%4.4i%2.2i.nc',pthdat,en4v,YR,IM);

    fprintf('EN4: Reading %i/%i \n',YR,IM);

    S = double(nc_varget(fnm,'salinity'));
    T = double(nc_varget(fnm,'temperature'))-273.15; % K -> C
    cc=cc+1;

    S = squeeze(S(1,:,:,:));
    T = squeeze(T(1,:,:,:));
    [ll,mm,nn]=size(S);

    clear SS TT
    if i2<i1,
      S1=S(:,j1:j2,i1:end);
      S2=S(:,j1:j2,1:i2);
      n1=size(S1,3);
      n2=size(S2,3);
      SS=S1;
      for in=1:n2
	n1=n1+1;
	SS(:,:,n1)=S2(:,:,in);
      end
      T1=T(:,j1:j2,i1:end);
      T2=T(:,j1:j2,1:i2);
      n1=size(T1,3);
      n2=size(T2,3);
      TT=T1;
      for in=1:n2
	n1=n1+1;
	TT(:,:,n1)=T2(:,:,in);
      end
    end

    sgm=sw_dens0(SS,TT)-1000;
    [ls,ms,ns]=size(sgm);
    mns=ms*ns;
    sgm=reshape(sgm,[ls,mns]);
    clear Zs0
    for im=1:mns
      s0=sgm(:,im);
      s0=[s0(1)-0.0001;s0];
      in=min(find(isnan(s0))); % interpolate below the bottom
      if ~isempty(in)
	s0(in:end)=s0(in-1);
      end

      z0=[0;ZZ];
      Zi=[0:-1:-1000]';
      s0i=interp1(z0,s0,Zi);
      iz=max(find(s0i<=sgm0));
      if ~isempty(iz),
	Zs0(im,1)=Zi(iz);
      else
	Zs0(im,1)=0;
      end
    end

    SGMD.TM(cc,1)=datenum(YR,IM,1);
    SGMD.AreaMean_Depth(cc,1) = nanmean(nanmean(Zs0));
    SGMD.Shallow_Depth(cc,1)  = max(max(Zs0));
    fprintf('Shallowest depth = %6.2fm\n',max(max(Zs0)));

  end % Time

  if s_mat>0
    fprintf('Saving %s\n',fmat);
    save(fmat,'SGMD');
  end
else
  fprintf('Loading %s\n',fmat);
  load(fmat);
end

TM = SGMD.TM;
DV = datevec(TM);
Dm = SGMD.AreaMean_Depth;
Dsh= SGMD.Shallow_Depth;

figure(1);
Im=find(DV(:,2)>=10 | DV(:,2)<=3); % winter
Is=find(DV(:,2)>=4 & DV(:,2)<=9); % summer


dw=Dm(Im);
ds=Dm(Is);

dw(dw<-400)=nan;
bind=[-400:20:0];
[Nw,Xw]=hist(dw,bind);

ds(ds<-400)=nan;
bind=[-400:20:0];
[Ns,Xs]=hist(ds,bind);

figure(2); clf;
axes('Position',[0.08 0.53 0.4 0.35]);
bar(Xw,Nw);
set(gca,'xlim',[-410 1],...
	'tickdir','out',...
	'ylim',[0 27],...
	'xtick',bind);
ttl=sprintf('GreenSea, AreaMn Wint(Oct-March), Z sgm=%4.2f, EN4, %i-%i',...
	    sgm0,yr1,yr2);
title(ttl);

axes('Position',[0.08 0.08 0.4 0.35]);
bar(Xs,Ns);
set(gca,'xlim',[-410 1],...
	'tickdir','out',...
	'ylim',[0 27],...
	'xtick',bind);
ttl=sprintf('GreenSea, AreaMn Sum(Apr-Sept), Z sgm=%4.2f, EN4, %i-%i',...
	    sgm0,yr1,yr2);
title(ttl);
set(gca,'xlim',[-410 1]);


% Boxplot
% Mean depth clim:
[a1,a2]=size(Dm);
nyr=a1/12;
Dd=reshape(Dm,[12,nyr]);
Dcl=mean(Dd,2);  

%Sdv=std(Dd')';  

axes('Position',[0.55 0.53 0.4 0.35]);
hold on;
%plot(Dcl,'k.-','linewidth',2);
%for ik=1:12
%  d1=Dcl(ik)-Sdv(ik);
%  d2=Dcl(ik)+Sdv(ik);
%  plot([ik ik],[d1 d2],'r-','Color',[0.6 0.6 0.6],'linewidth',2);
%end
hb=boxplot(Dd');
set(hb,'LineWidth',2);
ho=findobj(gca,'tag','Outliers');
set(ho,'MarkerEdgeColor',[1 1 1]);

set(gca,'ylim',[-600 -120],...
	'tickdir','out',...
	'xlim',[0.2 12.9],...
	'xtick',[1:12],...
	'ytick',[-550:50:100]);

ttl=sprintf('Gr. Sea, Clim., depths sgm=28, EN4, %i-%i',...
	    yr1,yr2);
title(ttl);

bottom_text(txtb);

  
if s_fig>0
  fgnm=sprintf('%sdepth_sgm%4.4i_EN4_GrSea',pthfig,round(sgm0*100));
  fprintf('Saving %s\n',fgnm);
  print('-depsc2',fgnm);
end


