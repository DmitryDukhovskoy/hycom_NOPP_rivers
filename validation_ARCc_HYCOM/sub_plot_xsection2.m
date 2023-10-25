function sub_plot_xsection2(nf,xx,zz,A,title_str,fld0,f_layer,...
		     xl1,xl2,yl1,yl2,c1,c2);
%
ZZ=zz;
plr=20; % highlight this interface in HYCOM

% Colormap:
nint= 200;
CMP = create_colormap5(nint,c1,c2);
cmp0= CMP.colormap;
for k=1:10
  cmp0(k,:)=[1 1 1];
end
nav = 15;
cmp = smooth_colormap(cmp0,nav);
%cmp=cmp0;
cnt = CMP.intervals;


figure(nf);
clf
% 
axes('position',[0.1 0.1 0.8 0.8]);
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
dmm=zz(k,I);
zb=min(dmm);
fprintf('Bottom zb=%6.1f\n',zb);

pcolor(xx,ZZ,A); shading interp;
%pcolor(xx,zz,A); shading flat;
caxis([c1 c2]);
hold on;
colormap(cmp);

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

x=xx(1,:);
set(gca,'Color',[0 0 0],...
	'tickdir','out',...
        'xlim',[xl1 xl2],...
        'ylim',[yl1 yl2]);
%set(gca,'xtick',[0:1000:max(x)],...
%	'ytick',[-4500:500:0]);
set(gca,'fontsize',12);


%  lt0=mean(LAT(j1,i1:i2));
title(title_str);
%text(xx(8),200,title_str,'Fontsize',16);

hb = colorbar;
set(hb,'TickLength',0.03);



return