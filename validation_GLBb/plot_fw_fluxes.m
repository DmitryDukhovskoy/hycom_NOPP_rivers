% Plot Calculate FW fluxes - annual mean
% across straits 
% /nexsan/GLBa0.08/expt_90.9/data
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

%s_mat  = 0; % 
s_fig  = 1;
f_pltS =1; % plot S in cross-sections


rg=9806;  % convert pressure to depth, m
%Sref=34.8;
Sref=35;

pth190='/nexsan/GLBb0.08/GLBb0.08_190/data/meanstd/'; % 1993-1994
pth191='/nexsan/GLBb0.08/GLBb0.08_191/data/meanstd/'; % 1995-2012
pthmat = '/Net/mars/ddmitry/hycom/GLBb0.08/data_mat/';
ftopo = sprintf('%sGLBb_T07_subset_Natl.mat',pthmat);
pthglb = '/nexsan/GLBb0.08/GLBb0.08_191/topo/';
pthfig = '/Net/mars/ddmitry/hycom/GLBb0.08/fig_fw/';
%fmat = sprintf('%sFWflx_annual_1993-2015.mat',pthmat);
fmat = sprintf('%sFWflx_annual_Sref%3i_1993-2015.mat',pthmat,round(Sref*10));
txtb='ddmitry/.../hycom_NOPP_rivers/validation_GLBb/plot_fw_fluxes.m';

fprintf('Loading FW fluxes %s\n',fmat);
load(fmat);
fprintf('Loading topo %s\n',ftopo);
load(ftopo);

[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2

nyrs=length(FWFLX);
for iy=1:nyrs
  YR(iy,1)=FWFLX(iy).Year;
end

Sref=FWFLX(1).Sref;

nsgm=7;
for iyr=1:nyrs
  for isgm=1:nsgm
    fwf=FWFLX(iyr).SEGM(isgm).FWFlux_m3;
    FWF(isgm,iyr)=nansum(nansum(fwf));
%    keyboard;
  end
end

% Plot sections
f_sect=0;
if f_sect>0
  cmp0=[0 0 0; 1 1 1];
  LMSK=HH*0;
  LMSK(HH<0)=1;
  [mm,nn]=size(HH);

  figure(10); clf;
  pcolor(LMSK); shading flat;
  colormap(cmp0);
  hold on;
  freezeColors;

  contour(HH,[-5000:1000:-50],'Linewidth',1.5,'Color',[0.7 0.7 0.7]);
  hold on
  for isgm=1:nsgm
    i1=FWFLX(1).SEGM(isgm).IJ(1,1);
    i2=FWFLX(1).SEGM(isgm).IJ(2,1);
    j1=FWFLX(1).SEGM(isgm).IJ(1,2);
    j2=FWFLX(1).SEGM(isgm).IJ(2,2);

    plot([i1 i2],[j1 j2],'r','Linewidth',2);
    ttl=sprintf('%i',isgm);
    text(0.5*(i1+i2),j1+10,ttl,'Fontsize',14);
  end
  axis('equal');
  set(gca,'xlim',[1 nn],...
	  'ylim',[1 mm],...
	  'xtick',[],...
	  'ytick',[]);

  nf=10;
  dlmb=20;
  dphi=10;
  clr=[0.8 0.8 0.8];
  plot_gridlines(dlmb,dphi, nf, clr, LON, LAT)

  bottom_text(txtb);
end

FWF(7,2)=nan; % some wierd jump ? in N.Atl section in 1994

POS(1,:)=[0.1 0.7 0.8 0.22];
POS(2,:)=[0.1 0.4 0.8 0.22];
POS(3,:)=[0.1 0.08 0.8 0.22];
%POS(4,:)=[0.08 0.05 0.8 0.15];
np=0;
ifg=1;
figure(ifg);clf;
for isgm=1:7
  np=np+1;
  if np>3
    bottom_text(txtb);
    ifg=ifg+1;
    figure(ifg); clf;
    np=1;
  end
  
  ps=POS(np,:);
  axes('Position',ps);
%  fwf=FWF(isgm,:);  % m3/s
  fwf=FWF(isgm,:)*3600*24*365*1e-9;  % km3/yr
%  yl1=floor(min(fwf));
%  yl2=ceil(max(fwf));
%  plot(YR,fwf,'Linewidth',2);
%  set(gca,'xlim',[1992.5 2015.5],...
%	  'xtick',[1993:2:2015],...
%	  'xminortick','on',...
%	  'tickdir','out',...
%	  'Fontsize',14);
  hb=bar(YR,fwf,0.9);
  set(hb,'Facecolor',[0.3 0.6 0.8]);
  yl1=floor(min(fwf)-0.1*abs(min(fwf)));
  yl2=ceil(max(fwf)+0.1*max(fwf));
  yl2=max([10,yl2]);
  yl1=min([-10,yl1]);
  set(gca,'xlim',[1995.5 2015.5],...
	  'xtick',[1993:1:2015],...
	  'ylim',[yl1 yl2],...
	  'ytick',[-10000:2000:10000],...
	  'tickdir','out',...
	  'fontsize',14);
  if isgm==6
    set(gca,'ytick',[-2000:200:600]);
  end
  
  ttl=sprintf('Annual mean FWFlux, km3/yr, Sref=%4.1f Section %i',Sref,isgm);
  title(ttl);
end

if s_fig>0
  for ii=1:ifg
    fgnm=sprintf('%sFWtransport_sections%2.2i',pthfig,ii);
    fprintf('Saving %s\n',fgnm);
    kll=sprintf('-f%i',ii);
    print('-depsc2',kll,fgnm);
  end
end

% Plot S fields in the cross sections:
% specify years to plot and section #
YRPLT=[1995; 2015];
isgm=1;
cyr=0;
clear Ssum Zsum
if f_pltS>0
% Long-term mean:
% Need to interpolate into fixed depths
  for iyr=1:nyrs
    ss=FWFLX(iyr).SEGM(isgm).S_section;
    zm=FWFLX(iyr).SEGM(isgm).ZM_section;
    IJ=FWFLX(iyr).SEGM(isgm).IJ;
    X=LON(IJ(1,2),IJ(1,1):IJ(2,1));
    hb=HH(IJ(1,2),IJ(1,1):IJ(2,1));
    year=YR(iyr);
%    ss=[ss(1,:);ss];
%    zm=[zm(1,:);zm];
%    zm(1,:)=0; % for interpolation
    [nl,npb]=size(ss);
    
    ZZi=[(0:-5:-500),(-510:-10:-1100),(-1150:-50:-2000),(-2100:-100:-3500)]';
    nli=length(ZZi);
    for ipp=1:npb
      Ib=min(find(zm(:,ipp)==0));
      if Ib==1, 
	Si(1:nli,ipp)=nan;
	continue; 
      end; % Land point
      s0=ss(:,ipp);
      z0=zm(:,ipp);
      dz=abs(diff(z0));
      ibz = min(find(dz<0.01))+1; % bottom
      if hb(ipp)<z0(ibz-1);
	z0(ibz-1)=hb(ipp); % match the bottom
      end
	
      for kl=ibz:nl
	z0(kl)=z0(kl-1)-0.1;
      end;
      z0=[0;z0];
      s0=[s0(1);s0];
      if z0(nl)>ZZi(end)
	z0(nl)=ZZi(end);
      end
      ibZi=max(find(ZZi>=hb(ipp)));
      s0i = interp1(z0,s0,ZZi,'cubic');
      s0i(ibZi:end)=nan;
      Si(:,ipp)=s0i;
    end

    
    if ~exist('Ssum')
      Ssum=Si*0;
      Zsum=zm*0;
    end
    Ssum=Ssum+Si;
    dmm=YRPLT-year;
    if ~all(dmm) 
      cyr=cyr+1;
      SS(cyr,:,:)=Si;
    end
    
    
  end; % year
  Smn = Ssum./nyrs;

  
% Plot S profiles:
% all settings are for Baffin bay section:
  nint=200;
  c1=33;
  c2=35;
  CMP = create_colormap2(nint,c1,c2);
  cnt=CMP.intervals;
  cmp=CMP.colormap;

  nm=FWFLX(1).SEGM(isgm).Name;
  ttl=sprintf('S %i %s',YRPLT(1),nm);

  figure(5); clf;
  S1=squeeze(SS(1,:,:));
  axes('position',[0.08 0.6 0.35 0.3]);
  pcolor(X,ZZi,S1); shading interp;
  colormap(cmp);
  title(ttl,'Fontsize',14);
  set(gca,'xlim',[-76.2 -64],...
	  'ylim',[-1500 0],...
	  'tickdir','out',...
	  'ytick',[-2000:250:0],...
	  'xtick',[-76:2:-62],...
	  'Fontsize',14,...
	  'Color',[0 0 0]);
  
  S1=squeeze(SS(2,:,:));
  ttl=sprintf('S %i %s',YRPLT(2),nm);
  axes('position',[0.5 0.6 0.35 0.3]);
  pcolor(X,ZZi,S1); shading interp;
  colormap(cmp);
  title(ttl,'Fontsize',14);
  set(gca,'xlim',[-76.2 -64],...
	  'ylim',[-1500 0],...
	  'tickdir','out',...
	  'ytick',[-2000:250:0],...
	  'xtick',[-76:2:-62],...
	  'Fontsize',14,...
	  'Color',[0 0 0]);
  

  hght=[];
  lngth=[];
  mint=20;
  mbx=mint;
  fsz=13;
  bxc='k';
  posc=[0.9 0.11 0.8 0.06];
  aend=0;
  [az,axc]  = colorbar_vert (cmp,cnt,hght,lngth,mint,fsz,bxc,posc,mbx,aend);
  
    bottom_text(txtb);
    
if s_fig>0
    fgnm=sprintf('%sS_section_%s',pthfig,nm);
    fprintf('Saving %s\n',fgnm);
%    kll=sprintf('-f%i',ii);
    print('-dpng','-f5','-r250',fgnm);
end

  
end





