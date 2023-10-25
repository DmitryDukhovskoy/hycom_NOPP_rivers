function sub_plot_fwc(YR,FWC,pthfig,s_fig,HH,LON,LAT);
% YR - array of year to plot
%      NxM, if M=1 - no long-term averaging
%           if M=2, for each n=1:N, FWC is averaged
%             for years yr1- yr2
fprintf('Plotting FWC ...\n');
[a1,a2]=size(YR);
if a2>2
  error('YR array can have 1 or 2 years in each row');
end

c1  = 0;
c2  = 30;
nint= 200;
CMP = create_colormap2_1(nint,c1,c2);
cmp0= CMP.colormap;
nav = 20;
cmp = smooth_colormap(cmp0,nav);
cnt = CMP.intervals;

avg=logical(0);
if a2==2
  avg=logical(1);
end

Nyr=length(FWC);
for k=1:Nyr
  YY(k,1)=FWC(k).Year;
end

[mm,nn]=size(HH);
Lm=HH;
Lm(Lm<0)=nan;
Lm(Lm>=0)=1;
cbb=[0 0 0; 1 1 1];

for k=1:a1  % years to plot
  yr1=YR(k,1);
  i1=find(YY==yr1);
  yr2=yr1;
  if avg
    yr2=YR(k,2);

    if yr2<yr1
      error('Year 2 (%i) has to be > Year 1 (%i)',yr2, yr1);
    end
    i2=find(YY==yr2);
    Fwc=zeros(mm,nn);
    cc=0;
    for ik=i1:i2
      cc=cc+1;
      Fwc=Fwc+FWC(ik).Fwc_m;
    end
    Fwc=Fwc./cc;
  else
    yr1=YR(k,1);
    Fwc=FWC(k).Fwc_m;
  end
  
  nf=1;
  figure(nf); clf;
  pcolor(Lm); shading flat;
  colormap(cbb);
  hold on;
  freezeColors;
  
  pcolor(Fwc); shading flat;
  caxis([c1 c2]);
  colormap(cmp);
  
%  contour(HH,[-200 -200],'Linewidth',1,...
%	  'Color',[0.9 0.9 0.9]);
  
  axis('equal');
  set(gca,'Color',[0 0 0]);
  set(gca,'xlim',[250 nn],...
	'ylim',[400 1950],...
	'xtick',[],...
	'ytick',[]);
  
  dlmb=20;
  dphi=10;
  clr=[0.8 0.8 0.8];
  plot_gridlines(dlmb,dphi, nf, clr, LON, LAT)
  
  hght=[];
  lngth=[];
  mint=20;
  mbx=mint;
  fsz=13;
  bxc='k';
  posc=[0.8 0.11 0.8 0.08];
  aend=0;
  [az,axc]  = colorbar_vert (cmp,cnt,hght,lngth,mint,fsz,bxc,posc,mbx,aend);

  ctl=sprintf('GLBb0.08, FWC (m), Sref=34.8, %i-%i',yr1,yr2);
  title(ctl,'Fontsize',13);
  
  if s_fig>0
    fgnm=sprintf('%sGLBb008_FWC_%i-%i',pthfig,yr1,yr2);
    fprintf('Saving %s\n',fgnm);
    print('-dpng','-r200',fgnm);
  end

  
end  % years to plot  
  
  



return