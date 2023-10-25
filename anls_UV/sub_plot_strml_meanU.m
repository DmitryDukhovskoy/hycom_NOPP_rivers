function sub_plot_strml_meanU(expt,yrF,arrow_dst,smx_L,smn_L,...
			      dlim,dip,strm_clr,v_col,lhead);
% Plot streamlines prepared
% in streamline_meanUV.m

%regn = 'ARCc0.08';
%expt = 110;
%yrF  = 2005

pthmat = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
fmat = sprintf('%sNAtlGreenl_particles_%i.mat',pthmat,yrF);

fprintf('Loading %s\n',fmat);
load(fmat);

fprintf('Plotting streamlines ...\n');

nd = length(PRTCL.TRACK);
for id=1:nd
  I = PRTCL.TRACK(id).I';
  J = PRTCL.TRACK(id).J';
  IP(id,:)=I;
  JP(id,:)=J;
end

np=size(IP,2);


%hold on;
%contour(HH,[-5000:1000:-10],'Color',[0.8 0.8 0.8]);

cpp=0;
cll=0;
XL=[]; % strml
YL=[];
XP=[]; % arrowhead
YP=[];
%arrow_dst = 1.5;
%dlim = 0.25; % min dist between isolines
for ip=1:dip:np
  
  if mod(cll,50)==0 & cll>0,
    fprintf('...   %6.3f%%\n',ip/np*100);
  end	
    
  
  X0a=IP(:,ip);
  Y0a=JP(:,ip);
  
  X0=X0a(~isnan(X0a));
  Y0=Y0a(~isnan(Y0a));
  dx=diff(X0);
  dy=diff(Y0);
  dsgm=sqrt(dx.^2+dy.^2);
  X0=X0(dsgm>0.1);
  Y0=Y0(dsgm>0.1);
  if length(X0)<round(1/10*length(X0a)), continue; end
  dx=diff(X0);
  dy=diff(Y0);
  dsgm=sqrt(dx.^2+dy.^2);
  dst=cumsum(dsgm);

  if dst(end)<smn_L, continue; end;
% Total length should be < smx_L  
  imx = max(find(dst<=smx_L));
  if imx<4, imx=4; end;
  if ~isempty(imx);
    X0=X0(1:imx);
    Y0=Y0(1:imx);
    dst2=dst(1:imx);
  end
% Find mean pnt:
  imean = max(find(dst2<=dst2(end)/2));
  if imean<3, imean=3; end;
  x1m = X0(imean-1);
  y1m = Y0(imean-1);
  x2m = X0(imean);
  y2m = Y0(imean);

% Do parameteric spline interpolation
  [Xi,Yi]=sub_parametric_spline(X0,Y0,100);
  Xi=Xi(:);
  Yi=Yi(:);

  X0=Xi;
  Y0=Yi;
%JJs=round(Yi);
%IIs=round(Xi);
  
  
  ncc=length(X0);
%  nc2=round(ncc/2);
%  dx1=((X0(1)-X0(end)).^2+(Y0(1)-Y0(end)).^2);
%  dx2=((X0(1)-X0(nc2)).^2+(Y0(1)-Y0(nc2)).^2);
%  if dx1<0.1 & dx2<0.1, continue; end; % short streamline
% Check if the streamlines are too crowded:
  if isempty(XL), 
    XL=X0'*1e9;
    YL=Y0'*1e9;
  end
  nx=length(Xi);
  nmn=round(nx/2);
  xmn=X0(nmn);
  ymn=Y0(nmn);
% end  
  DL = ((X0(end)-XL).^2+(Y0(end)-YL).^2);
  dlmn = min(min(DL));
  if dlmn<dlim, continue; end;
% start
  DL = ((X0(1)-XL).^2+(Y0(1)-YL).^2);
  dlmn = min(min(DL));
  if dlmn<dlim, continue; end;
% middle
  DL = ((xmn-XL).^2+(ymn-YL).^2);
  dlmn = min(min(DL));
  if dlmn<dlim, continue; end;
  cll=cll+1;
  XL(cll,:)=X0';
  YL(cll,:)=Y0';
  
  plot(X0,Y0,'-','Color',strm_clr);
% Plot arowhead  - scaled to be same size
  x2 = x2m;
  y2 = y2m;
  if lhead==0
    x1=x1m;
    y1=y1m;
  else 
    a1=(x2m-x1m);
    a2=(y2m-y1m);
    aa=sqrt(a1^2+a2^2);
    a1m = a1/aa*lhead;
    a2m = a2/aa*lhead;
    x1 = x2m-a1m;
    y1 = y2m-a2m;
  end
  
  if isempty(XP), XP=X0'*1e9; YP=Y0'*1e9; end;
  D=((X0(end)-XP).^2+(Y0(end)-YP).^2);
  dmn=min(D);
  if dmn>arrow_dst
    cpp=cpp+1;
    XP(cpp,:)=X0';
    YP(cpp,:)=Y0';
    cf=1;
    beta=15;
%    v_col = [0 0.4 0.6];
    lwd=1;
    draw_arrowF(x1,x2,y1,y2,cf,beta,v_col,lwd);
%    keyboard
  end
  
end

return


