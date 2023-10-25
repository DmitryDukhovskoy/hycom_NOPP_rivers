function sub_plot_xsct(nfg,ZZ,XL,F,pfld,dy,sname);

switch(pfld);
 case('temp')
  c1=-1.8;
  c2=4.8;
  dc = 0.4;
  nint = 200;
%  CMP = create_colormap5(nint,c1,c2);
%  cmp = CMP.colormap;
%  cnt = CMP.intervals;
%
%  nav = 15;
%  cmp = smooth_colormap(cmp,nav);
  cmp = colormap(jet(200));
 case('saln');
  c1=33;
  c2=35;
  dc=0.25;
%  cmp = colormap(parula(200));
  cmp = colormap(jet(200)); % to match with Woodgate et al, 2005
  
end



figure(nfg); clf;
AX=axes('Position',[0.11 0.12 0.77 0.8]);
pcolor(XL,ZZ,F); shading interp
hold on
%pcolor(XL,ZZ,F); shading flat
switch(sname)
 case('BarOp');
  if c1>30
    c1=34.5;
    c2=35.1;
    dc=0.1;
  else
    c1=-1;
    c2=8;
    dc=0.5;
  end
  set(gca,'xdir','reverse');
 case('BeringS')
  if c1>30
%    c1=31.4;
%    c2=32.4;
%    dc=0.1;
    c1=28.75;
    c2=33.25;
    dc=0.25;
  else
%    c1=-1;
%    c2=4;

    c1=1.5;
    c2=11.5;
    dc=0.5;
  end
%  set(gca,'xdir','reverse');
 case('FramS');
  if c1>20
    c1=28.5;
    c2=35.5;
    dc =1;
  else
    c1=-1.8;
    c2=4.8;
    dc=0.4;
  end
 case('DavisS');
  if c1>30
    c1=32.5;
    c2=34.5;
    dc=0.25;
  else
    c1=-2;
    c2=2;
    dc=0.2;
  end
    
end


caxis([-12000 12000]);
caxis([c1 c2]);
colormap(AX,cmp);
hc = colorbar;
set(hc,'Position',[0.89 0.11 0.025 0.8],...
       'Ticks',[c1:dc:c2],...
       'TickLength',0.042,...
       'Fontsize',11);
%keyboard
%stt=sprintf('%s, Volume Flux, m3/s',nm);
%title(stt);
xl1=XL(1);
xl2=XL(end);
if xl1>xl2
  xl1=xl2;
  xl2=XL(1);
end

set(gca,'Color',[0 0 0],...
	'tickdir','out',...
	'xlim',[xl1 xl2],...
	'ylim',[min(min(ZZ))-1 0],...
	'ytick',[-5000:dy:0]);


return