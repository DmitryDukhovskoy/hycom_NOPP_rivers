  function sub_plot_fld(A,fgn,ttl,HH);
% Quick plot of sea ice field

figure(fgn); clf;
pcolor(A); shading flat;
hold on;

if (~isempty(HH));
  contour(HH,[0 0],'k');
  caxis([nanmin(nanmin(A)) nanmax(nanmax(A))]);
end;

colorbar;

title(ttl);



return
