function UV = sub_set_xsct_v2(HH,LON,LAT,SEGM);
% Define vertical sections
% specified in SEGM
% sections are either along X or Y 
nsgm = length(SEGM);
for ig=1:nsgm
  segm=SEGM{ig};
  fprintf('%i, %s\n',ig,segm);
  switch(segm)
   case('BeringS');
    is1=634;
    js1=1919;
    is2=658;
    js2=js1;
    X=LON(js1,is1:is2);
    Y=LAT(js1,is1:is2);
    I=[is1:is2]';
    J=ones(size(I))*js1;
   case('FramS');
%    is1=935; too far north - ~80N
%    js1=959;
%    is2=1070;
%    js2=js1;
    is1=908; % this ~78.5
    is2=1106;
    js1=915;
    js2=915;
    X=LON(js1,is1:is2);
    Y=LAT(js1,is1:is2);
    I=[is1:is2]';
    J=ones(size(I))*js1;
   case('BarOp');
    is1=1166;
    js1=456;
    is2=is1;
    js2=929;
    X=LON(js1:js2,is1);
    Y=LAT(js1:js2,is1);
    J=[js1:js2]';
    I=ones(size(J))*is1;
   case('DavisS');
    is1=478;
    js1=680;
    is2=568;
    js2=680;
    X=LON(js1,is1:is2);
    Y=LAT(js1,is1:is2);
    I=[is1:is2]';
    J=ones(size(I))*js1;
   case('DenmarkS');
    is1=705;
    js1=550;
    is2=855;
    js2=js1;
    X=LON(js1,is1:is2);
    Y=LAT(js1,is1:is2);
    I=[is1:is2]';
    J=ones(size(I))*js1;
   case('IclNorw');
    is1=899;
    js1=477;
    is2=1176;
    js2=js1;
    X=LON(js1,is1:is2);
    Y=LAT(js1,is1:is2);
    I=[is1:is2]';
    J=ones(size(I))*js1;
   case('FaroeShetl');
    is1=1032;
    js1=1099;
    is2=429;
    js2=js1;
    X=LON(js1,is1:is2);
    Y=LAT(js1,is1:is2);
    I=[is1:is2]';
    J=ones(size(I))*js1;
   case('NorthEGr');
    is1=908;
    js1=843;
    is2=1004;
    js2=js1;
    if is1>is2, a=is2; is2=is1; is1=a; end;
    X=LON(js1,is1:is2);
    Y=LAT(js1,is1:is2);
    I=[is1:is2]';
    J=ones(size(I))*js1;
   case('CntrEGr');
    is1=874;
    js1=668;
    is2=900;
    js2=js1;
    if is1>is2, a=is2; is2=is1; is1=a; end;
    X=LON(js1,is1:is2);
    Y=LAT(js1,is1:is2);
    I=[is1:is2]';
    J=ones(size(I))*js1;
   case('SouthEGr');
    is1=628;
    js1=463;
    is2=650;
    js2=js1;
    if is1>is2, a=is2; is2=is1; is1=a; end;
    X=LON(js1,is1:is2);
    Y=LAT(js1,is1:is2);
    I=[is1:is2]';
    J=ones(size(I))*js1;
   case('SouthWGr');
    is1=531;
    js1=500;
    is2=554;
    js2=js1;
    if is1>is2, a=is2; is2=is1; is1=a; end;
    X=LON(js1,is1:is2);
    Y=LAT(js1,is1:is2);
    I=[is1:is2]';
    J=ones(size(I))*js1;
   case('CntrWGr');
    is1=553;
    js1=802;
    is2=605;
    js2=js1;
    if is1>is2, a=is2; is2=is1; is1=a; end;
    X=LON(js1,is1:is2);
    Y=LAT(js1,is1:is2);
    I=[is1:is2]';
    J=ones(size(I))*js1;
   case('NarresS');
    is1=628;
    js1=1034;
    is2=is1;
    js2=1054;
%    if is1>is2, a=is2; is2=is1; is1=a; end;
    X=LON(js1:js2,is1);
    Y=LAT(js1:js2,is1);
    J=[js1:js2]';
    I=ones(size(J))*is1;
   case('LancasterS');
    is1=475;
    js1=1018;
    is2=507;
    js2=js1;
    if is1>is2, a=is2; is2=is1; is1=a; end;
    X=LON(js1,is1:is2);
    Y=LAT(js1,is1:is2);
    I=[is1:is2]';
    J=ones(size(I))*js1;
  end
  
  UV(ig).Name = segm;
  UV(ig).X    = X;
  UV(ig).Y    = Y;
  UV(ig).I    = I;
  UV(ig).J    = J;
  ll=length(I);
  clear D
  
  for k=1:ll-1
    x1=X(k);
    y1=Y(k);
    x2=X(k+1);
    y2=Y(k+1);
    D(k,1)=distance_spheric_coord(y1,x1,y2,x2);
  end
%keyboard
  D(k+1,1)=D(k);
  UV(ig).Dist = D;

end

  p_map=0;
  if p_map>0
%    ig = 2;
    figure(11); clf;
    contour(HH,[0 0],'k');
    hold on
    contour(HH,[-1000 -1000],'c');
    for ig=1:nsgm
      I = UV(ig).I;
      J = UV(ig).J;
      plot(I,J,'r.-');
    end
    keyboard
  end



return