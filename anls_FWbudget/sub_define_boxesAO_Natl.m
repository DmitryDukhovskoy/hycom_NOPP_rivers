function BX = sub_define_boxesAO_Natl(HH,LON,LAT,fplot);
%
% Define regions ("boxes") in the interior
% N.Atl convective sites and Arctic Ocean
% based on box coordinates
% This function can be used on any
% LON/LAT grid as long as specified boxes
% are inside the domain
%

fprintf('Finding indices for study regions ...\n');

[mm,nn]=size(HH);
[II,JJ] = meshgrid((1:nn),(1:mm));

cc=0;
% Arctic Ocean use getpts
XY = [         -170.745483398438          66.3834838867188
         -166.311218261719          65.6099319458008
           -157.3134765625           66.302978515625
         -133.643432617188          68.1148452758789
          -111.91357421875          66.3755874633789
         -94.6777954101562                 66.171875
            -87.1455078125          66.7610931396484
         -82.1927185058594          73.0902557373047
         -81.2518615722656          75.0712509155273
         -68.8072814941406          77.9717025756836
         -19.5343933105469          79.9050903320312
          13.1079711914062          79.2272033691406
          22.8864135742188          69.2546234130859
          36.4021606445312          61.5729522705078
          68.9031982421875          67.5582656860352
               82.99609375           67.892463684082
          126.119232177734           69.539680480957
         -177.362670898438          66.8302307128906];

IJ = sub_XY2indx(XY,LON,LAT);
cc=cc+1;
BX(cc).Name = 'ArcticOcean';
BX(cc).XY = XY;
BX(cc).IJ = IJ;
IN = inpolygon(II,JJ,IJ(:,1),IJ(:,2));
BX(cc).IN_polygon=find(IN==1);

% interior Labrador sea
XY = [-54.84, 56.48; ...
      -49.01, 57.11; ...
      -50.38, 60.33; ...
      -56.33, 59.49];
IJ = sub_XY2indx(XY,LON,LAT);
cc=cc+1;
BX(cc).Name = 'LabradorSea';
BX(cc).XY = XY;
BX(cc).IJ = IJ;
IN = inpolygon(II,JJ,IJ(:,1),IJ(:,2));
BX(cc).IN_polygon=find(IN==1);

% Irminger
XY = [-43.87, 56.85;...
      -36.27, 56.55; ...
      -33.96, 60.44; ...
      -40.12, 60.44];
IJ = sub_XY2indx(XY,LON,LAT);
cc=cc+1;
BX(cc).Name = 'IrmingerSea';
BX(cc).XY = XY;
BX(cc).IJ = IJ;
IN = inpolygon(II,JJ,IJ(:,1),IJ(:,2));
BX(cc).IN_polygon=find(IN==1);

% Iceland
XY = [-16.44, 68.16; ...
      -9.06, 68.07; ...
      -8.43, 70.71; ...
      -16.49, 70.895];
IJ = sub_XY2indx(XY,LON,LAT);
cc=cc+1;
BX(cc).Name = 'IcelandSea';
BX(cc).XY = XY;
BX(cc).IJ = IJ;
IN = inpolygon(II,JJ,IJ(:,1),IJ(:,2));
BX(cc).IN_polygon=find(IN==1);

% Greenland
XY = [-5.41, 73; ...
      2.33, 72.57; ...
      4.812, 75.104; ...
      -3.91, 75.626];
IJ = sub_XY2indx(XY,LON,LAT);
cc=cc+1;
BX(cc).Name = 'GreenlandSea';
BX(cc).XY = XY;
BX(cc).IJ = IJ;
IN = inpolygon(II,JJ,IJ(:,1),IJ(:,2));
BX(cc).IN_polygon=find(IN==1);

% Eurasian Basin
XY =[-37.65, 86.28; ...
     137.15, 79.44; ...
     113.95, 78.68; ...
     12.46, 82.29];
IJ = sub_XY2indx(XY,LON,LAT);
cc=cc+1;
BX(cc).Name = 'EurasianBasin';
BX(cc).XY = XY;
BX(cc).IJ = IJ;
IN = inpolygon(II,JJ,IJ(:,1),IJ(:,2));
BX(cc).IN_polygon=find(IN==1);

% Canadian Basin - BeaufSea
XY =[-132.31, 73.09; ...
     -157.49, 72.91; ...
     -159.07, 83.24; ...
     -116.84, 81.75];
IJ = sub_XY2indx(XY,LON,LAT);
cc=cc+1;
BX(cc).Name = 'BeaufortSea';
BX(cc).XY = XY;
BX(cc).IJ = IJ;
IN = inpolygon(II,JJ,IJ(:,1),IJ(:,2));
BX(cc).IN_polygon=find(IN==1);

if fplot>0
  fprintf('Plotting Domain with %i regions\n',cc);
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