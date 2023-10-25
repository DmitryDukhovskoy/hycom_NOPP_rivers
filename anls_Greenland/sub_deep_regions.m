function BX = sub_deep_regions(HH,LON,LAT,h0,fplot);
%
% Define deep basins/regions in the interior
% Subpolar N. Atlantic
%
%[a,b]=getpts
%IJ=round([a,b]) 
%I=sub2ind(size(HH),IJ(:,2),IJ(:,1));
%[LON(I),LAT(I)]

fprintf('Finding indices for deep regions ...\n');
cc=0;
% interior Labrador sea
XY = [-68.8934020996094          66.7906951904297
      -63.2389526367188           65.761116027832
         -51.8563232421875          67.3300476074219
         -43.6701049804688          60.3960418701172
         -57.3036499023438          52.2764930725098
          -65.410888671875          58.4321784973145];
IJ = sub_XY2indx(XY,LON,LAT);
cc=cc+1;
BX(cc).Name = 'LabradorSea';
BX(cc).XY = XY;
BX(cc).IJ = IJ;
fprintf('%i, %s\n',cc,BX(cc).Name);

% Irminger
XY = [-43.6701049804688          60.3960418701172
      -42.2982788085938          65.1557998657227
      -33.931396484375          68.2893447875977
      -23.8016357421875          65.5860290527344
      -21.3866271972656          65.2572860717773
      -35.4088134765625          55.2817916870117];
IJ = sub_XY2indx(XY,LON,LAT);
cc=cc+1;
BX(cc).Name = 'IrmingerSea';
BX(cc).XY = XY;
BX(cc).IJ = IJ;
fprintf('%i, %s\n',cc,BX(cc).Name);

% Iceland
XY = [-33.931396484375          68.2893447875977
      -26.0924987792969          73.5082244873047
      -8.348388671875          71.0728530883789
      -7.14614868164062          62.0742492675781
      -21.3866271972656          65.2572860717773
      -23.8016357421875          65.5860290527344];
IJ = sub_XY2indx(XY,LON,LAT);
cc=cc+1;
BX(cc).Name = 'IcelandSea';
BX(cc).XY = XY;
BX(cc).IJ = IJ;
fprintf('%i, %s\n',cc,BX(cc).Name);

% Greenland
XY = [-26.0924987792969          73.5082244873047
      -16.4093017578125          81.1646270751953
       12.3865356445312          79.3284301757812
       19.0155639648438          74.4336929321289
      -8.348388671875          71.0728530883789];
IJ = sub_XY2indx(XY,LON,LAT);
cc=cc+1;
BX(cc).Name = 'GreenlandSea';
BX(cc).XY = XY;
BX(cc).IJ = IJ;
fprintf('%i, %s\n',cc,BX(cc).Name);

% Norwegian
XY = [19.0155639648438          74.4336929321289
      21.1515197753906          69.2700958251953
      6.43115234375          59.7063980102539
     -7.14614868164062          62.0742492675781
     -8.348388671875          71.0728530883789];
IJ = sub_XY2indx(XY,LON,LAT);
cc=cc+1;
BX(cc).Name = 'NorwegianSea';
BX(cc).XY = XY;
BX(cc).IJ = IJ;
fprintf('%i, %s\n',cc,BX(cc).Name);

%[a,b]=getpts
%IJ=round([a,b]) 
%I=sub2ind(size(HH),IJ(:,2),IJ(:,1));
%[LON(I),LAT(I)]

% Baffin Bay
XY =[-63.2389526367188           65.761116027832
     -70.9468994140625          68.7379455566406
     -82.0050659179688          71.7297515869141
     -82.405517578125          73.2074279785156
     -81.0233459472656          77.1445999145508
     -63.8454284667969           78.242561340332
     -49.4643249511719          70.8137664794922
     -51.8563232421875          67.3300476074219];
IJ = sub_XY2indx(XY,LON,LAT);
cc=cc+1;
BX(cc).Name = 'BaffinBay';
BX(cc).XY = XY;
BX(cc).IJ = IJ;
fprintf('%i, %s\n',cc,BX(cc).Name);

%IJ=round([a,b]) 
%I=sub2ind(size(HH),IJ(:,2),IJ(:,1));
%[LON(I),LAT(I)]
% Central and Eastern Subpolar Gyre
XY=[-57.3036499023438          52.2764930725098
    -43.6701049804688          60.3960418701172
    -35.4088134765625          55.2817916870117
    -21.3866271972656          65.2572860717773
    -7.14614868164062          62.0742492675781
    6.43115234375          59.7063980102539
       6.83642578125          50.7333526611328
   -1.27902221679688          48.1440734863281
   -56.3408508300781          48.1440734863281];
IJ = sub_XY2indx(XY,LON,LAT);
cc=cc+1;
BX(cc).Name = 'CntrEast Subpolar Gyre';
BX(cc).XY = XY;
BX(cc).IJ = IJ;
fprintf('%i, %s\n',cc,BX(cc).Name);


% Select deep domain within each region
%h0=-800;
[mm,nn]=size(HH);
[XX,YY] = meshgrid((1:nn),(1:mm));
nbx=cc;
for ib=1:nbx
  iBG = BX(ib).IJ;
  INP = inpolygon(XX,YY,iBG(:,1),iBG(:,2));
  IN = find(INP==1 & HH<h0);
  BX(ib).IN = IN;
%  fprintf(' Tracer integrated for region: %i %s\n',ib,BX(ib).Name);
end
  

CLR=[0.8 0.8 1;
     0.8 1 1;
     1 0.9 0.9;
     1 1 0.8;
     1 0.9 1;
     0.9 0.95 0.9;
     0.95 1 0.9;
     0.85 1 0.9;
     0.95 0.7 1];

if fplot>0
  fprintf('Plotting Domain with %i regions\n',cc);
  figure(10); clf;
  contour(HH,[0 0],'k');
  hold on;
  contour(HH,[-200 -200],'Color',[0.8 0.8 0.8]);
  contour(HH,[-2000 -2000],'Color',[0.8 0.8 0.8]);

  cc = size(BX,2);
  for ik=1:cc
    IN=BX(ik).IN;
    clr=CLR(ik,:);
    plot(XX(IN),YY(IN),'.','Color',clr);
    
    IJ=BX(ik).IJ;
    IJ(end+1,:)=IJ(1,:);
    plot(IJ(:,1),IJ(:,2),'r.-');
    
    x0=mean(IJ(:,1));
    y0=mean(IJ(:,2));
    spp=sprintf('# %i',ik);
    text(x0,y0,spp,'Fontsize',14);
  end
  
  axis('equal');
  set(gca,'tickdir','out',...
	  'xlim',[320 1300],...
	  'ylim',[100 1100]);
  title('SPNA deep regions FWC analysis');
  tbt = 'sub_deep_regions.m';
  bottom_text(tbt,'pwd',1);
end

  
return