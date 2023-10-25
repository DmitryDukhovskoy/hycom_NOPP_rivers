% HYCOM 0.04:
%
% Analyze characteristics of the 
% coastal currents on the Southeastern
% Greenland shelf EGCC and EGC
% following
% Sutherland and Pickart, PiO, 2008
% Lentz and Largier, JPO, 2005
% 
% SE Gr. Shelf vertical transsections
% extracted in flux_SEGrShelf_xsct004.m
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

f_mat=1;

%Sc0=33.5; % salinity defining plume front
Sc0=33.8; % salinity defining plume front
%Sc0=33.6; % salinity defining plume front
%Sc0=34.; % salinity defining plume front
%Sc0=33.7; % salinity defining plume front

regn = 'ARCc0.04';
%expt = 010;
expt = 012; % Greenland runoff
ys=2005;
ye=2006;

pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/topo_grid/';
pthmat = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/%3.3i/data_GrSect/',expt);

ftopo = sprintf('%s/depth_%s_17DD.nc',pthtopo,regn); % 
btx='anls_EGC004.m';

fprintf('%s-%3.3i, Salin Front=%4.2f\n',regn,expt,Sc0);

% If need map of section - see flux_SHGrShelf_xsct.m

% Select figures to plot:
f_dltS=1;  % bottom - surf S
f_FWC=1;   % FWC
f_RW=1;    % Ratio of onshore/offshore plume transect areas
f_brich=0; % Bulk Rich.
f_Sbtanm=1;
nfg=0;

fprintf('Analyzing EGC front, Salinity=%4.2f\n',Sc0);
fprintf('The following figs turned on: \n');
if f_dltS==1
  fprintf('     Bottom S - Surf S\n');
end
if f_FWC==1
  fprintf('     FWC on the shelf\n');
end
if f_RW==1
  fprintf('     Ratio of plume areas\n');
end
if f_brich==1
  fprintf('     Bulk Rich Number\n');
end



f_plt=0;
if f_plt==1
  nfg=1;
  dnmb=datenum(2000,08,4);
  dv=datevec(dnmb);
  YR=dv(1);
  fmat=sprintf('%sarc08_expt%3.3i_SEGreenlSh_xsct_dayUTSZ_%i.mat',pthmat,expt,YR);
  xl1=0;
  xl2=90;
  yl1=-1000;
  yl2=0;
  
  sub_plot_GrShelf_xsect(dnmb,nfg,fmat,xl1,xl2,yl1,yl2,nm,Sc0);
  bottom_text(btx,'pwd',1);
end

% Find distance of the front foot and offshore distance on surface
% and calculate plume cross-sectional areas 
% Aon - onshore from the max depth (foot)
% Aoff - offshore
foutp=sprintf('%sarc04_%3.3i_anls_EGC%3.3i.mat',pthmat,expt,Sc0*10);

cc=0;
for iyr=ys:ye
  if f_mat==2; break; end;
  YR=iyr;
%  fmat=sprintf('%sarc04_expt%3.3i_SEGreenlSh_xsct_dayUTSZ_%i.mat',pthmat,expt,YR);
  fmat=sprintf('%s%3.3i_SEGreenlSh_xsct_dayUTSZ_%i.mat',pthmat,expt,YR);
  fprintf('Loading %s\n',fmat);
  load(fmat);
  
  TM=UTSZ.Time;
  ZZ=UTSZ.ZZlevels;
  DZ=abs(diff(ZZ));
  ZM=UTSZ.ZM;
  nzm=length(ZM);
  nzz=length(ZZ);
  LL=UTSZ.Dist_origin; % m
  DL=UTSZ.Segm_dL; % m - segment lengths
  nm=UTSZ.Name;
  Hb=UTSZ.Hbottom;
  Hb(Hb>=0)=0;

  nrec=length(TM);
%
% Get indices of the bottom cells
  if ~exist('Ibtm','var')
    S=squeeze(UTSZ.Saln(1,:,:));
    [mm,nn]=size(S);
    ii1=min(find(~isnan(S(1,:))));
    Ibtm=S(1,:)*nan;
    for kk=ii1:nn
      k2=max(find(~isnan(S(:,kk))));
      if isempty(k2); continue; end;
      Ibtm(kk)=sub2ind(size(S),k2,kk);
    end
%    Ibtm=Ibtm(~isnan(Ibtm));
  end
 
% Calculate mean S:
  S=UTSZ.Saln;
  Smean=squeeze(nanmean(S,1));
  inn=find(~isnan(Ibtm));
  Sbt_mean=Smean(1,:)*0;
  Sbt_mean(inn)=Smean(Ibtm(inn));
  
  
  
  for ik=1:nrec
    S=squeeze(UTSZ.Saln(ik,:,:));
    [mm,nn]=size(S);
  
% Get bottom S:
    inn=find(~isnan(Ibtm));
    Sb=S(1,:)*nan;
    Sb(inn)=S(Ibtm(inn));
% To avoid nans in searching min S, make 1st land point Smin
    Sb(1)=min([Sc0-0.05, Sb(2)-0.01]);
    
    Ic=min(find(Sb>=Sc0));
    if isempty(Ic),
      fprintf('Cannot locate front %4.2s\n',Sc0);
      return
    end
    
    dc=LL(Ic-1)+(LL(Ic)-LL(Ic-1))/(Sb(Ic)-Sb(Ic-1))*(Sc0-Sb(Ic-1));
% Max depth of the front - where it intersects the
% bottom, foot of the front:
    Hp=Hb(Ic-1)+(Hb(Ic)-Hb(Ic-1))/(Sb(Ic)-Sb(Ic-1))*(Sc0-Sb(Ic-1));
% Area vrt xsection, onshore part:
    Aon=abs(nansum(Hb(1:Ic-1).*DL(1:Ic-1)))+...
	abs((dc-LL(Ic-1))*Hp);
    
% surface
    Ss=S(1,:);
    Is=min(find(Ss>=Sc0));
    ds=LL(Is-1)+(LL(Is)-LL(Is-1))/(Ss(Is)-Ss(Is-1))*(Sc0-Ss(Is-1));
    Is=Is-1;
% Area vrt xsection, offshore part:
% find depth of Sc0
    Aoff=0;
    if Ic==Is
      Aoff=abs(Hb(Is)*(ds-dc));
    elseif Is<Ic | ds<dc
      Aoff=1e-1;
    else
     for kk=Ic+1:Is
      k1=max(find(S(:,kk)<=Sc0));
      if isempty(k1), break; end;
      if isnan(S(k1+1,kk))  % bottom
	zz0=Hb(kk);
      else
        zz0=ZM(k1)+(ZM(k1+1)-ZM(k1))/(S(k1+1,kk)-S(k1,kk))*(Sc0-S(k1,kk));
      end
      
      Aoff=Aoff+abs(zz0*DL(kk));
     end
    end
% Difference: btm-surf S:
    dS=Sb-Ss;

% Bottom S anomaly:
    Sanom=Sb-Sbt_mean;

%
% FWC:
    Sref=34.8;
    dltS=(Sref-S)./Sref;
    dltS(isnan(dltS))=0;
    FWC=(dltS'*DZ)./abs(Hb)*50; % m of FW in 50 m of water

% Bulk Richardson Number (approximation of grad. Rich.#)
% <0.25 - turbilent regime
    T=squeeze(UTSZ.Temp(ik,:,:));
    U=squeeze(UTSZ.Unrm(ik,:,:));
    Rho=sw_dens0(S,T);
    Rib=zeros(1,nn)*nan;
    Umn=zeros(1,nn)*nan;
    for iss=1:nn
      if isnan(S(1,iss)) | S(1,iss)>Sc0, continue; end;
      kS=max(find(S(:,iss)<=Sc0));
      rho1=nanmean(Rho(1:kS,iss));
      Hsal=abs(sum(DZ(1:kS)));
      u1=nansum(U(1:kS,iss).*DZ(1:kS))/Hsal; % depth-averaged U

      Umn(iss)=u1;
      
      if isnan(S(kS+1,iss)) % bottom
	kS=kS-1;
      end
      
      rho2=Rho(kS+1,iss);
      dltR=rho2-rho1;
      u2=U(kS+1,iss);
      dltU=(u2-u1);
      rho0=1028;
      Rib(iss)=9.8*dltR.*abs(ZZ(kS))./(rho0*dltU.^2);
    end
    U_mean=nanmean(Umn);
    
    cc=cc+1;
    PLM.S_front=Sc0;
    PLM.TM(cc)=TM(ik);
    PLM.DistBtm(cc)=dc;
    PLM.DistSrf(cc)=ds;
    PLM.dltSvrt(cc,:)=dS;
    PLM.Sbtm_anom(cc,:)=Sanom;
    PLM.Aonshore(cc)=Aon;
    PLM.Aoffshore(cc)=Aoff;
    PLM.FWC_m50m(cc,:)=FWC';
    PLM.BulkRich(cc,:)=Rib;
    PLM.Front_foot_depth(cc)=Hp;
    PLM.Umean_plume(cc)=U_mean;
    PLM.Hbottom=Hb;
    
    
  end % days 
  
  fprintf('%i finished\n',iyr);
  
  if f_mat==1
    fprintf('Saving %s\n',foutp);
    save(foutp,'PLM');
  end
  
  
end

if f_mat==2
  fprintf('Loading %s\n',foutp);
  load(foutp);
  
end

% Plotting:

TM=PLM.TM;
Aon=PLM.Aonshore;
Aoff=PLM.Aoffshore;
Rw=Aon./Aoff;
%Rw(Rw>1)=1;
dS=PLM.dltSvrt;
FWC=PLM.FWC_m50m;
Sbanm=PLM.Sbtm_anom;

DV=datevec(TM);
ndy=length(find(DV(:,1)==DV(1,1)));
YRS=[DV(1,1):1/ndy:DV(end,1)+0.9999];

ctt=0;
ytck=[];
for yr=DV(1,1):DV(end,1)
  for im=1:12
    ii=find(DV(:,1)==yr & DV(:,2)==im,1);
    ctt=ctt+1;
    ytck(ctt)=YRS(ii);
    ylbl{ctt}=sprintf('%2.2i/%4.4i',im,yr);
  end
end


% Plot bottom-surf S difference
LM=LL*1e-3;

if f_dltS==1
  cl1=colormap_red(100);
  cl2=flipud(colormap_blue(100));
  cmp=[cl2;cl1];
  c1=-0.6;
  c2=0.6;

  nfg=nfg+1;
  figure(nfg); clf;
  axes('Position',[0.2 0.1 0.25 0.8]);
  pcolor(LM,YRS,dS); shading interp;
  colormap(cmp);
  caxis([c1 c2]);
  hold on;
  plot(PLM.DistBtm*1e-3,YRS,'-','Color',[0 0 0],'linewidth',2);
  plot(PLM.DistSrf*1e-3,YRS,'-','Color',[1 1 0],'linewidth',2);

  cb=colorbar;
  set(cb,'Position',[0.47 0.25 0.02 0.5],...
	 'Ticks',[c1:0.1:c2],...
	 'TickLength',0.04,...
	 'Fontsize',14);
  xl1=0;
  xl2=60;
  set(gca,'tickdir','out',...
	  'xlim',[0 60],...
	  'xtick',[0:10:60],...
	  'ytick',ytck,...
	  'yticklabel',ylbl,...
	  'Fontsize',14);
  xlabel('Distance offshore, km');

  stl=sprintf('%s-%3.3i, %s, Sbtm-Ssrf',regn,expt,nm);
  title(stl,'Interpreter','none');

  axes('Position',[0.5 0.1 0.1 0.05]);
  hold on;
  plot([0.1 0.3],[0.1 0.1],'k-','Linewidth',2);
  plot([0.1 0.3],[0.2 0.2],'k-','Color',[1 1 0],'Linewidth',2);
  text(0.35,0.1,sprintf('%4.2f bottom',Sc0),'Fontsize',14);
  text(0.35,0.2,sprintf('%4.2f surf',Sc0),'Fontsize',14);
  set(gca,'Visible','off',...
	  'xlim',[0 1],...
	  'ylim',[0.08 0.23]);

  % Plot bottom profile along transect
  axes('Position',[0.62 0.78 0.25 0.1])
  hbx=[LM(1);LM;LM(end)];
  hby=[-5000;Hb;-5000];
  fill(hbx,hby,[0 0 0]);

  set(gca,'tickdir','out',...
	  'xlim',[xl1 xl2],...
	  'ylim',[-500 0],...
	  'xtick',[0:10:xl2],...
	  'ytick',[-4000:100:0],...
	  'Fontsize',12);
  xlabel('Distance, km');

  bottom_text(btx,'pwd',1);

  set(gcf,'Position',[932 39 1259 1292]);

end


% Plot Area onshore/Area offshore
% of the plume
cc1=0;
cc2=2;
dmm=colormap_WB(200,cc1,cc2);
cmp2=dmm.colormap;

if f_FWC==1
  nfg=nfg+1;
  figure(nfg); clf;
  axes('Position',[0.2 0.1 0.25 0.8]);
  pcolor(LM,YRS,FWC); shading interp;
  colormap(cmp2);
  caxis([cc1 cc2]);
  hold on;
  plot(PLM.DistBtm*1e-3,YRS,'-','Color',[0 0 0],'linewidth',2);
  plot(PLM.DistSrf*1e-3,YRS,'-','Color',[1 1 0],'linewidth',2);


  cb=colorbar;
  set(cb,'Position',[0.47 0.25 0.02 0.5],...
	 'Ticks',[0:0.2:2],...
	 'TickLength',0.04,...
	 'Fontsize',14);

  set(gca,'tickdir','out',...
	  'xlim',[0 60],...
	  'xtick',[0:10:60],...
	  'ytick',ytck,...
	  'yticklabel',ylbl,...
	  'Fontsize',14);
  xlabel('Distance offshore, km');

  stl=sprintf('arc004-%3.3i, %s, FWC (m in 50m), Sref=34.8',expt,nm);
  title(stl,'Interpreter','none');
  xlabel('Distance offshore, km');

  axes('Position',[0.5 0.1 0.1 0.05]);
  hold on;
  plot([0.1 0.3],[0.1 0.1],'k-','Linewidth',2);
  plot([0.1 0.3],[0.2 0.2],'k-','Color',[1 1 0],'Linewidth',2);
  text(0.35,0.1,sprintf('%4.2f bottom',Sc0),'Fontsize',14);
  text(0.35,0.2,sprintf('%4.2f surf',Sc0),'Fontsize',14);
  set(gca,'Visible','off',...
	  'xlim',[0 1],...
	  'ylim',[0.08 0.23]);

  % Plot bottom profile along transect
  axes('Position',[0.62 0.78 0.25 0.1])
  hbx=[LM(1);LM;LM(end)];
  hby=[-5000;Hb;-5000];
  fill(hbx,hby,[0 0 0]);

  set(gca,'tickdir','out',...
	  'xlim',[xl1 xl2],...
	  'ylim',[-500 0],...
	  'xtick',[0:10:xl2],...
	  'ytick',[-4000:100:0],...
	  'Fontsize',12);
  xlabel('Distance, km');


  bottom_text(btx,'pwd',1);
  %axes('Position',[

  set(gcf,'Position',[930 39 1259 1292]);

end

% ---------------------------
%  PLOT Bulk Rich Number
% ---------------------------
Rib = PLM.BulkRich;
if f_brich==1
  cc1=0;
  cc2=2;
  dmm=colormap_WB(200,cc1,cc2);
  cmp2=dmm.colormap;

  nfg=nfg+1;
  figure(nfg); clf;
  axes('Position',[0.2 0.1 0.25 0.8]);
  pcolor(LM,YRS,Rib); shading interp;
  colormap(cmp2);
  caxis([cc1 cc2]);
  hold on;
  plot(PLM.DistBtm*1e-3,YRS,'-','Color',[0 0 0],'linewidth',2);
  plot(PLM.DistSrf*1e-3,YRS,'-','Color',[1 1 0],'linewidth',2);

  set(gca,'tickdir','out',...
	  'xlim',[0 60],...
	  'xtick',[0:10:60],...
	  'ytick',ytck,...
	  'yticklabel',ylbl,...
	  'Fontsize',14);
  xlabel('Distance offshore, km');

end

%  ===============================
%  Plotting Area onshore/offshore
%  ===============================
if f_RW==1
  nfg=nfg+1;
  
  figure(nfg); clf;
  axes('Position',[0.2 0.1 0.25 0.8]);
  plot(log10(Rw),YRS);
  hold on;
  plot([0 0],[YRS(1) YRS(end)],'k--');

  set(gca,'tickdir','out',...
	'xlim',[-2 8],...
	'xtick',[-2:8],...
	'ylim',[YRS(1) YRS(end)],...
	'ytick',ytck,...
	'yticklabel',ylbl,...
	'xgrid','on',...
	'ygrid','on',...
	'Fontsize',14);
  stl=sprintf('arc004-%3.3i, %s, log(Ab/As)',expt,nm);
  title(stl,'Interpreter','none');
  xlabel('log10(Ab/As)');

  axes('Position',[0.5 0.8 0.3 0.1]);
  clear stl
  stl{1}=sprintf('Ratio Ab/As=onshelf/offshelf area FW current');
  stl{2}=sprintf('log10 scale, when Ab/As<1 is ');
  stl{3}=sprintf('buoyant current surface trapped');
  stl{4}=sprintf('Front is defined by S=%4.2f\n',Sc0);
  text(0,0,stl,'Fontsize',12);
  set(gca,'Visible','off');

  bottom_text(btx,'pwd',1);
  %axes('Position',[

%  set(gcf,'Position',[930 39 1259 1292]);

end


%  ===============
%  S bttom anomaly
%  =============
if f_Sbtanm==1
  cb1=-0.3;
  cb2=0.3;
  dmm=colormap_WB(200,cb1,cb2);
  a1=flipud(dmm.colormap);
  dmm=colormap_WG(200,cb1,cb2);
  a2=dmm.colormap;
  cmp3=[a1;a2];
  
  
  nfg=nfg+1;
  figure(nfg); clf;
  axes('Position',[0.2 0.1 0.25 0.8]);
  pcolor(LM,YRS,Sbanm); shading interp;
  colormap(cmp3);
  caxis([cb1 cb2]);
  hold on;
  plot(PLM.DistBtm*1e-3,YRS,'-','Color',[0 0 0],'linewidth',2);
  plot(PLM.DistSrf*1e-3,YRS,'-','Color',[0.7 0 0],'linewidth',2);


  cb=colorbar;
  set(cb,'Position',[0.47 0.25 0.02 0.5],...
	 'Ticks',[cb1:0.1:cb2],...
	 'TickLength',0.04,...
	 'Fontsize',14);

  set(gca,'tickdir','out',...
	  'xlim',[0 60],...
	  'xtick',[0:10:60],...
	  'ytick',ytck,...
	  'yticklabel',ylbl,...
	  'Fontsize',14);
  xlabel('Distance offshore, km');

  stl=sprintf('arc004-%3.3i, %s, Sbtm anomaly',expt,nm);
  title(stl,'Interpreter','none');
  xlabel('Distance offshore, km');

  axes('Position',[0.5 0.1 0.1 0.05]);
  hold on;
  plot([0.1 0.3],[0.1 0.1],'k-','Linewidth',2);
  plot([0.1 0.3],[0.2 0.2],'k-','Color',[0.7 0 0],'Linewidth',2);
  text(0.35,0.1,sprintf('%4.2f bottom',Sc0),'Fontsize',14);
  text(0.35,0.2,sprintf('%4.2f surf',Sc0),'Fontsize',14);
  set(gca,'Visible','off',...
	  'xlim',[0 1],...
	  'ylim',[0.08 0.23]);

  % Plot bottom profile along transect
  axes('Position',[0.62 0.78 0.25 0.1])
  hbx=[LM(1);LM;LM(end)];
  hby=[-5000;Hb;-5000];
  fill(hbx,hby,[0 0 0]);

  set(gca,'tickdir','out',...
	  'xlim',[xl1 xl2],...
	  'ylim',[-500 0],...
	  'xtick',[0:10:xl2],...
	  'ytick',[-4000:100:0],...
	  'Fontsize',12);
  xlabel('Distance, km');


  bottom_text(btx,'pwd',1);
end


