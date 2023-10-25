function BX = sub_define_boxes(HH,LON,LAT,fplot);
%
% Define regions ("boxes") in the interior
% N.Atl convective sites and Arctic Ocean
% based on box coordinates
% This function can be used on any
% LON/LAT grid as long as specified boxes
% are inside the domain
%

fprintf('Finding indices for study regions ...\n');
cc=0;
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

% Baffin Bay
XY =[-73.3792, 73.6356; ...
     -63.397, 68.962; ...
     -59.629, 69.842; ...
     -70.118, 74.658];
IJ = sub_XY2indx(XY,LON,LAT);
cc=cc+1;
BX(cc).Name = 'BaffinBay';
BX(cc).XY = XY;
BX(cc).IJ = IJ;

% Subpolar Gyre:
% Labr - Irm - eastern & central N. Atlantic
XY=[-67.395751953125          66.8995208740234
         -62.5514831542969          66.2832336425781
         -52.7092895507812          67.8802185058594
         -28.0538635253906          68.8377914428711
         -22.7550964355469          65.5771102905273
         -15.6484069824219          64.5237350463867
          -6.6981201171875          62.1011924743652
         -1.02557373046875          60.4161567687988
          6.84304809570312          59.6347846984863
          9.00018310546875          56.0720176696777
          7.43930053710938          50.0580825805664
         -3.67416381835938          48.2008666992188
         -54.5737915039062          47.9641075134277
         -57.2241516113281          51.8668746948242
         -64.7034606933594          59.4301719665527
         -66.7424926757812          62.2638473510742
         -71.3219604492188          65.7100067138672];
IJ = sub_XY2indx(XY,LON,LAT);
cc=cc+1;
BX(cc).Name = 'Subpolar Gyre';
BX(cc).XY = XY;
BX(cc).IJ = IJ;

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