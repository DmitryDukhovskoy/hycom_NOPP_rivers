   function sub_plot_divu(DIV,LON,LAT,nf,HH,U,V,stl,xl1,xl2,yl1,yl2,c1,c2);
%
%
%

% Colormap:
na=100;
cl1=colormap_gray(na);
cl2=colormap_purple(na);
cl3=colormap_blue(na);
cl4=colormap_green(na);
cl5=colormap_yellow(na);
cl6=colormap_red(na);

%cmp=[cl1;cl2;cl3;cl4;cl5;cl6];
%cmp=[flipud(cl3);cl6];
%CMP=colormap_blue_cyan_white(200,c1,c2);
%cmp=CMP.colormap;
%cnt=CMP.intervals;
%for ik=1:10
%  cmp(ik,:)=[1 1 1];
%end
%cmp=smooth_colormap(cmp,7);
cmp=colormap(parula(200));
nint=length(cmp);
%c1=-0.2;
%c2=0.2;
cnt=(c1:(c2-c1)/nint:c2);  % construct array of intervals

DIV(:,1:xl1)=nan;
DIV(:,xl2:end)=nan;
DIV(1:yl1,:)=nan;
DIV(yl2:end,:)=nan;
figure(nf);
clf
pcolor(DIV); shading flat;
hold on;
contour(HH,[0 0],'k','Linewidth',1);
contour(HH,[-5000:1000:-100],'Color',[0.7 0.7 0.7],'Linewidth',1);
%caxis([0 2]);
caxis([c1 c2]);
colormap(cmp);
axis('equal');
%  set(gca,'xlim',[0 nn],'ylim',[0 mm]);
set(gca,'xlim',[xl1 xl2],'ylim',[yl1 yl2]);
set(gca,'xtick',[],'ytick',[]);
clr=[0.9 0.9 0.9];
plot_gridlines(10,15,1,clr,LON,LAT);
%stl=sprintf('divU %4.4i/%2.2i/%2.2i, Layer %i',DV(1:3),plr);
title(stl,'Fontsize',12);

f_vct=1;
if isempty(U),
  f_vct=0;
end

if f_vct>0
  fprintf('Plotting vectors ...\n');
  scl=70;
  cf=0.5;
  beta=20;
  lwd=1.;
  v_col=[0 0 0];
  dii=round((xl2-xl1)/100);
  if dii==0, dii=1; end;
  nij=(xl2-xl1)/dii*(yl2-yl1)/dii;
  ncc=0;
  for ii=xl1:dii:xl2
    for jj=yl1:dii:yl2
      ncc=ncc+1;
      if (mod(ncc,500)==0)
	fprintf('... %4.1f%% done \n', ncc/nij*100);
      end
	    
      clear u v
      u = U(jj,ii);
      v = V(jj,ii);
      if isnan(u), continue; end;
      x0=ii;
      y0=jj;
      s0=sqrt(u^2+v^2);
      x1=x0+u*scl;
      y1=y0+v*scl;
      draw_arrowF(x0,x1,y0,y1,cf,beta,v_col,lwd);
    end
  end
end

%  keyboard


hght=[];
lngth=[];
mint=10;
mbx=mint;
fsz=12;
bxc='k';
posc=[0.85 0.1 0.8 0.06];
aend=1;
[az,axc]  = colorbar_vert(cmp,cnt,hght,lngth,mint,fsz,bxc,posc,mbx,aend);


%  freezeColors;
%  pcolor(hmsk); shading flat;
%  colormap([0 0 0]);

%colorbar
drawnow

return
  
