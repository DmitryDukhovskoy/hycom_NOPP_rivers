function BX = sub_define_Yash_reg(HH,LON,LAT,fplot);
% Prepare regions for analyzing 
% dlt S
% USing Yashayaev's regions
POBS = sub_obs_locations;
nobs = length(POBS);
for ii=1:nobs
  BX(ii).Name = POBS(ii).Name;
  X=POBS(ii).Lon;
  Y=POBS(ii).Lat;
  X=X(:);
  Y=Y(:);
  if length(X)==1
    x=X;
    y=Y;
    X=[x-0.16;...
       x+0.16;...
       x+0.16;...
       x-0.16;...
       x-0.16];
    Y=[y+0.16;...
       y+0.16;...
       y-0.16;...
       y-0.16;...
       y+0.16];
  end
  
  IJ = sub_indices_domain([X,Y],LON,LAT); 
  BX(ii).XY=[X,Y];
  BX(ii).IJ=IJ;
end  

%keyboard

if fplot>0
  fprintf('Plotting Domain with %i regions\n',ii);
  figure(10); clf;
  contour(HH,[0 0],'k');
  hold on;
  contour(HH,[-200 -200],'Color',[0.8 0.8 0.8]);
  contour(HH,[-2000 -2000],'Color',[0.8 0.8 0.8]);

  cc = size(BX,2);
  for ik=1:cc
    IJ=BX(ik).IJ;
    IJ(end+1,:)=IJ(1,:);
    plot(IJ(:,1),IJ(:,2),'r.-');
    x0=mean(IJ(:,1));
    y0=mean(IJ(:,2));
    spp=sprintf('# %i',ik);
    text(x0,y0,spp,'Fontsize',14);
    axis('equal');
    tbt = 'sub_define_boxes.m';
    bottom_text(tbt,'pwd',1);
  end
end



return