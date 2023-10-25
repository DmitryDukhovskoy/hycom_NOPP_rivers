function sub_plot_tracers(lTr,nf,HH,xlim1,xlim2,...
			  ylim1,ylim2,LON,LAT,stl,varargin);
%
% Can specify limits for plotting: 'c1',-6,'c2',3
%
if (nf<0),
  nf=abs(nf);
  figure('Visible','off'); clf;
  fprintf('Figure window is turned off\n');
else  
  figure(nf);
  clf
end


hmsk=HH;
hmsk(HH<0)=nan;

c1=[];
c2=[];
f_cmp = [];
nV = length(varargin);
%keyboard
for k=1:nV
  aa=varargin{k};
  aa=lower(aa);
  if strncmp(aa,'c1',2);
    c1=varargin{k+1};
  elseif strncmp(aa,'c2',2);
    c2=varargin{k+1};
  elseif strncmp(aa,'cmp',3);
    f_cmp=varargin{k+1};
  end;  
end

fprintf('Plotting %s ...\n',stl);

% Colormap:
%f_cmp=2;
if isempty(f_cmp),
  f_cmp=3;
end
if isempty(c1) | isempty(c2)
  c1=-3;
  c2=2;
%c1=-2;
%c2=4;
%if nTr==5 % Bering Str.
%  c1=-3;
%  c2=3;
%end
%if nTr==1 % Greenland
%  c1=-3;
%  c2=3;
%end
end
fprintf('    Limits climit: c1=%6.2f c2=%6.2f\n',c1,c2);
%keyboard
%c1=-2;
%c2=3;
switch(f_cmp)
 case(1)
  na=50;
%  cl1=colormap_gray(na);
  cl2=colormap_purple(na);
  cl2(1,:) = [1 1 1];
  cl2(2,:) = [1,1,1];
  cl2(3,:) = [1,1,1];
  cl3=colormap_blue(na);
  cl4=colormap_green(na);
  cl5=colormap_yellow(na);
  cl6=colormap_red(na);
  cmp=[cl2;cl3;cl4;cl5;cl6];
  cmp=smooth_colormap(cmp,10);
  nint=length(cmp);
  cnt=(c1:(c2-c1)/nint:c2);  % construct array of intervals
 case(2)
  nint= 240;
  CMP = create_colormap2_3(nint,c1,c2);
  cmp0= CMP.colormap;
  nav = 6;
  cmp = smooth_colormap(cmp0,nav);
  cnt = CMP.intervals;
 case(3)
  cmp = colormap(parula(200));
%  cmp(1,:) = [1 1 1];
  for k=1:20
    cmp(k,:)=[1 1 1];
  end
  cmp = smooth_colormap(cmp,15);
  cmp=flipud(cmp);
  nint= length(cmp);
  cnt=(c1:(c2-c1)/nint:c2);  % construct array of intervals
 case(4)
  cmp = colormap(parula(200));
  cmp(end,:)   = [1 1 1];
  cmp(end-1,:) = [1 1 1];
  cmp(end-2,:) = [1 1 1];
  cmp(end-3,:) = [1 1 1];
  nint= length(cmp);
  cmp = smooth_colormap(cmp,5);
  cmp(end,:)   = [1 1 1];
  cnt=(c1:(c2-c1)/nint:c2);  % construct array of intervals
 case(5)
  cmp = colormap_red(200);
  for k=1:10
    cmp(k,:)=[1 1 1];
  end
  cmp(end,:)=[0.5 0 0];
  cmp = smooth_colormap(cmp,5);
  nint= length(cmp);
  cnt=(c1:(c2-c1)/nint:c2);  % construct array of intervals
 case(6)
  cmp = colormap_blue(200);
  for k=1:20
    cmp(k,:)=[1 1 1];
  end
  cmp(end,:)=[0 0 0.4];
  cmp = smooth_colormap(cmp,15);
  nint= length(cmp);
  cnt=(c1:(c2-c1)/nint:c2);  % construct array of intervals
%  keyboard
 case(7)
  nint= 240;
  CMP = create_colormap_freshwater(nint,c1,c2);
  cmp0= CMP.colormap;
  nav = 6;
  cmp = smooth_colormap(cmp0,nav);
  cnt = CMP.intervals;
end

nint=length(cmp);


% Normalize:
%mT=max(max(Tr));
%Tr=Tr./mT*100;

pcolor(lTr); shading flat;
hold on;
%contour(HH,[0 0],'k','Linewidth',1);
contour(HH,[-5000:1000:-100],'Color',[0.7 0.7 0.7],'Linewidth',1);
%caxis([0 2]);
caxis([c1 c2]);
colormap(cmp);
axis('equal');
set(gca,'xlim',[xlim1 xlim2],'ylim',[ylim1 ylim2]);
set(gca,'xtick',[],'ytick',[]);
clr=[0.9 0.9 0.9];
plot_gridlines(45,10,1,clr,LON,LAT);
title(stl,'Fontsize',12,'Interpreter','none');

%hb = colorbar('EastOutside');
%set(hb,

hght=[];
lngth=[];
mint=20;
mbx=mint;
fsz=12;
bxc='k';
posc=[0.87 0.1 0.8 0.04];
aend=0;
[az,axc]  = colorbar_vert(cmp,cnt,hght,lngth,mint,fsz,bxc,posc,mbx,aend);

freezeColors;
pcolor(hmsk); shading flat;
colormap([0 0 0]);




return