% ARCc0.04
% Plot vertical disctribution
% of Greenland Tracer  in fjords
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/seawater;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

nTr = 1;
s_mat = 0; % =0 - extract data, no save; =1 extract & save, =2 - load saved
s_fig = 1;
s_map = 0; % plot map with segments
%f_tracloc = 0; % =1 - add locations of tracer release to transect maps
rg=9806;  % convert pressure to depth, m
%fld0='salin'; %'temp' - background field

plr=0; % highlight this interface
btx = 'xsection_tracer_fjords004.m';

regn = 'ARCc0.04';
expt = 011;
pthfig = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/%s/%3.3i/fig_green_xsct/',...
		  regn,expt);
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/topo_grid/';
%pthmat  = sprintf('/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
%fmat=sprintf('%svert_vel_layers.mat',pthmat);

ftopo = sprintf('%s/depth_%s_17DD.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);


YRPLT=[];
cc=0;
for iyr=2006:2006
  for idd=15:30:365
    cc=cc+1;
    YRPLT(cc,1)=iyr;
    YRPLT(cc,2)=idd;
  end
end

np   = size(YRPLT,1);
SGM  = sub_fjord_sections004(HH,LON,LAT);
nsgm = length(SGM);

f_map_check=0; % quickly see the segments
if f_map_check>0
  figure(10); clf;
%  fn = 10;
%  sub_plot_bath(HH,LON,LAT,fn,domname);
  hold on
  contour(HH,[0 0],'k');
  contour(HH,[-1000:100:-100],'Color',[0.8 0.8 0.8]);
  contour(HH,[-5000:1000:-1000],'Color',[0.9 0.9 0.9]);
  for kk=1:nsgm
    IIs = SGM(kk).IIs;
    JJs = SGM(kk).JJs;
    plot(IIs,JJs,'b.-');
    stx=sprintf('Sect %2.2i %s',kk,SGM(kk).Name);
    text(min(IIs),min(JJs),stx);
  end
  axis('equal');
  set(gca,'xlim',[1000 1920],...
	  'ylim',[800 2160]);
  bottom_text(btx,'pwd',1);
end

%fprintf('Section: %s, Saving fig: %i\n',xname,s_fig);
%keyboard

% Plot fields:
cnc=0;  % 
ip1=1;
for ip=ip1:np
  yr=YRPLT(ip,1);
  iday=YRPLT(ip,2);
%  pthbin = sprintf('/nexsan/archive/%s_%3.3i/data/%i/',regn,expt,yr)
  pthbin = sprintf('/nexsan/hycom/%s_%3.3i/data/%i/',regn,expt,yr)

  fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
  finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);

  
  if ~exist(fina,'file');
    fprintf('Not found: %s\n\n',fina);
    continue;
  end
  
  cnc=cnc+1;
  dnmb=datenum(yr,1,1)+iday-1;
  DV=datevec(dnmb);

  fprintf('%4.4i_%2.2i_%2.2i: %s\n',DV(1:3),fina);

  tic;
  [F,n,m,l] = read_hycom(fina,finb,'tracer','r_tracer',nTr);
  F(F>1e6)=nan;
  TR = F;

  fld='thknss';
  [F,n,m,l] = read_hycom(fina,finb,fld);
  F(F>1e18)=nan;
  F=F/rg;
  F(F<1e-2)=0;
  dH = F;
  
  for kk=1:nsgm
    INDs = SGM(kk).INDs;
    xnm  = SGM(kk).Name;
    Tr   = squeeze(TR(:,INDs));
  
% Prepare for plotting - add extra bogus layer
% at the bottom, Matlab won't show it
    [a1,a2]=size(Tr);
  
% Prepare vertical thicknesses for plotting  
    Dsec=squeeze(dH(:,INDs));
    Dsec(Dsec==0)=nan;
% Create Depth array of interface depths:
% Note these are BOTTOM interfaces 
% So Layer 1 is between interfaces 0m and ZZ(1)
    clear ZZb
    Dsec(isnan(Dsec))=0;
    ZZb(1,:)=-Dsec(1,:);
    for kl=2:l
      ZZb(kl,:)=ZZb(kl-1,:)-Dsec(kl,:);
    end
% For plotting need to have #of layers +1 interfaces
% add surface layer, otherwise all values
% will be shifted 1 layer down
    [nl,npb]=size(ZZb);
    ZZ=zeros(nl+1,npb);
    ZZ(2:nl+1,:)=ZZb;
    ZZav=ZZ;
  
% Depths of middle of the cells:
    clear ZM
    ZM(1,:)=0.5*ZZ(1,:);
    for kl=1:l
      ZM(kl,:)=0.5*(ZZ(kl+1,:)+ZZ(kl,:));
    end

% Interpolate into reg. vert. grid:
    ZZi=[0:-5:-1000]';
    nli=length(ZZi);
    Tri=zeros(nli,npb)*nan;
    Hs=HH(INDs);
    for ipp=1:npb
      Ib=min(find(Dsec(:,ipp)==0));
      if Ib==1, continue; end;
      t=Tr(:,ipp);
      zz=ZM(:,ipp);
      hb=-sum(Dsec(:,ipp));
      ibz = max(find(ZZi>=hb));
      for kl=Ib:nl
	zz(kl)=zz(kl-1)-0.1;
      end;
      zz=[0;zz];
      t=[t(1);t];
      if zz(nl)>ZZi(end)
	zz(nl)=ZZi(end);
      end
      ti = interp1(zz,t,ZZi,'pchip');
      ti(ibz+1:end) = nan;
      Tri(:,ipp) = ti;
    end

    XL = SGM(kk).dist_m*1e-3; % distance along segment, km  
    if strncmp(xnm,'WGr',3),
      xdr=-1;
    else
      xdr=1;
    end

% Plot:
    Tri(Tri<=0)=nan;
    lTr = log(Tri);
    
    fld0='tracer';
    stl=sprintf('%3.3i, trcr %i, sct %2.2i, %4.4i/%2.2i/%2.2i',...
	      expt, nTr, kk, DV(1:3));
    nf=1; % <0 visible off
    sub_plot_xsectZZ(nf,XL,ZZi,lTr,stl,fld0,xdr,Hs);
    
    if s_fig>0
      bottom_text(btx,'pwd',1);
      ffg=sprintf('%s%3.3i_sgm%2.2i_trcr%i_%4i%2.2i%2.2i',...
		pthfig,expt,kk,nTr,DV(1:3));
      fprintf('Saving %s\n\n',ffg);
      print('-dpng','-r200',ffg);
    end

  end % segments
  
  fprintf('1 time processed %6.1fmin\n\n',toc/60);
    
end



% --------------------
% Plot the cross-section
% delete unneeded topo:
if s_map==1
  HH(1:760,:)=nan;
  HH(2200:end,:)=nan;
  HH(:,2000:end)=nan;
  HH(:,1:800)=nan;

  fprintf('Plotting segments ...\n');
  fn=2;
  xyl = [930, 790;...
	 1900, 2200];
  
  sub_plot_bath2(HH,LON,LAT,fn,xyl);
  contour(HH,[-100:10:-5],'Color',[0.8 0.8 0.8]);
  contour(HH,[-3500:100:-10],'Color',[0.6 0.6 0.6]);
  contour(HH,[-4000:500:-50],'Color',[0.3 0.3 0.3]);
  for kk=1:nsgm
    IIs = SGM(kk).IIs;
    JJs = SGM(kk).JJs;
    x   = SGM(kk).dist_m*1e-3; % km
    plot(IIs,JJs,'b-','Linewidth',2);
    for km=1:20:max(x)
      d=abs(x-km);
      i0=find(d==min(d));
      plot(IIs(i0),JJs(i0),'r.','Markersize',9);
    end
    stx=sprintf('Sect %2.2i %s',kk,SGM(kk).Name);
    text(min(IIs),min(JJs),stx);
  end
  title('ARCc0.04 Fjord sections');

  bottom_text(btx,'pwd',1);


  tt = ' ';
  if f_tracloc>0
    ftrca = sprintf('%srelax_trcr_Greenland.a',PTH.data);
    ftrcb = sprintf('%srelax_trcr_Greenland.b',PTH.data);
    imonth = 7;
    FTr = sub_read_tracer_relax(ftrca,ftrcb,IND,imonth,nlev);
    lFTr = log(FTr); % kg/m3

    c1=-15;
    c2=5;
    nint=100;
    cnt=(c1:(c2-c1)/nint:c2);  % construct array of intervals
    Rr=-1;                        % colors of max intensity, clockwise rotation 
    C0=[-1,0,1];                 % starting point, red dimension is Z axis
    Cend=[1, 1, 0.6];
    cmpT = colormap_spiral(nint,'C0',C0,'Rr',Rr,'Cend',Cend);

    pcolor(lFTr); shading flat;
    colormap(cmpT);

    hght=[];
    lngth=[];
    mint=10;
    mbx=mint;
    fsz=12;
    bxc='k';
    posc=[0.85 0.1 0.8 0.06];
    aend=0;
    colorbar_vert(cmpT,cnt,hght,lngth,mint,fsz,bxc,posc,mbx,aend);
    tt=sprintf(', Tracer Rlx, kg/m3, mo=%i',imonth);

  end



  %axis('equal');
  %switch (xname)
  % case('uummq')
  %  set(gca,'xlim',[170 300],'ylim',[600 720]); 
  %  stl = sprintf('Fjord Uummannaq %s',tt);
  % case('scoresby');
  %  set(gca,'xlim',[440 580],'ylim',[500 630]); 
  %  stl = sprintf('Scoresby Sound%s',tt);
  %end  
  %set(gca,'xtick',[],'ytick',[]);
  %title(stl,'Fontsize',16);

  %bottom_text(txtb);

  ffg=sprintf('%sfjord_xsct_map',pthfig);
  if s_fig>0
    fprintf('Saving %s\n\n',ffg);
    print('-dpng','-r200',ffg);
  end

end
