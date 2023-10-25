function sub_plot_xsection3(nf,xx,zz,A,title_str,fld0,f_layer,...
		     xl1,xl2,yl1,yl2,c1,c2,cntr,s0);
% cntr - [] not contours
%

ZZ=zz;
plr=20; % highlight this interface in HYCOM

% Colormap:
nint= 200;
CMP = create_colormap5(nint,c1,c2);
cmp0= CMP.colormap;
for k=1:round(0.06*nint);
  cmp0(k,:)=[1 0.4 1];
end
nav = 15;
cmp = smooth_colormap(cmp0,nav);
%cmp=cmp0;
cnt = CMP.intervals;


figure(nf);
clf
% 
axes('position',[0.1 0.1 0.82 0.8]);
%xx=XL;
%zz=ZZav;
%A=lTr;
% Find bottom:
[ll,nn]=size(A);
for k=ll:-1:1
  sb=A(k,:);
  I=find(~isnan(sb));
  if ~isempty(I), break; end;
end
%keyboard
dmm=zz(k,I);
zb=min(dmm);
fprintf('Min Bottom zb=%6.1f\n',zb);

pcolor(xx,ZZ,A); shading interp;
%pcolor(xx,zz,A); shading flat;
caxis([c1 c2]);
hold on;
colormap(cmp);
%keyboard
if f_layer>0
  nl=size(zz,1);
  for k=1:nl
    z=ZZ(k,:);
    x=xx(k,:);
    plot(x,z,'k-','linewidth',1);
    if (k==plr+1), % bottom interface of the layer
      plot(x,z,'k-','linewidth',1.8);
    end
  end
end;

%keyboard
if ~isempty(cntr),
  contour(xx,ZZ,A,cntr,'k');
  contour(xx,ZZ,A,[cntr(1) cntr(1)],'k','Linewidth',1.6);
end


x=xx(1,:);
set(gca,'Color',[0 0 0],...
	'tickdir','out',...
        'xlim',[xl1 xl2],...
        'ylim',[yl1 yl2]);
%set(gca,'xtick',[0:1000:max(x)],...
%	'ytick',[-4500:500:0]);
set(gca,'fontsize',14);


%  lt0=mean(LAT(j1,i1:i2));
title(title_str);
%text(xx(8),200,title_str,'Fontsize',16);

hb = colorbar;

if strncmp(fld0,'salin',4)
  tkk=get(hb,'Ticks');
  for ik=1:length(tkk)
    stk(ik)=salin2exp(s0,tkk(ik),-1);
    stk(ik)=(round(stk(ik)*100))/100;
    lbstk{ik}=sprintf('%4.2f',stk(ik));
  end;

  set(hb,'TickLabels',lbstk);
%keyboard
end
set(hb,'TickLength',0.03,...
       'Position',[0.93 0.1 0.012 0.8],...
       'Fontsize',14);


set(gcf,'Position',[997 626 1530 705]);

return
