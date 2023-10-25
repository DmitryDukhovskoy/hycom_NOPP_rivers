function PRTCL = sub_seed_Greenl(HH,LON,LAT);
% Seed particles around Greenland
% in N. Atlantic
%
iBG = [         502        1033
         465        1001
         431         957
         409         719
         363         629
         377         582
         373         503
         379         267
         422         241
         406         160
         465         124
         620         113
         901         509
        1118         948
         962        1110
         723        1135
         665        1087
         557        1045
         545        1032
         515        1033];

[mm,nn]=size(LON);

[XX,YY] = meshgrid((1:nn),(1:mm));
INP = inpolygon(XX,YY,iBG(:,1),iBG(:,2));
XX(~INP)=nan;
YY(~INP)=nan;
XX(HH>-10)=nan;
YY(HH>-10)=nan;
dp=15;
XP=XX(1:dp:mm,1:dp:nn);
YP=YY(1:dp:mm,1:dp:nn);
I=find(~isnan(XP));
II=XP(I);
JJ=YP(I);

chck=0;
if chck>0
  figure(10); clf;
  contour(HH,[0 0],'k');
  hold on;
  contour(HH,[-1000 -1000],'b');
  plot(iBG(:,1),iBG(:,2),'r.-');
%  contour(LON,[-153 -153],'g');
%  contour(LON,[-130 -130],'g');
%  contour(LAT,[73 73],'g');
%  contour(LAT,[83 83],'g');
%  [JJ,II]=ind2sub(size(HH),IN);
  plot(II,JJ,'b.');
  axis('equal');
end


PRTCL.TRACK.I=[];
PRTCL.TRACK.J=[];

icnt=1;
nprt=0; % particle counter
  
PRTCL.TRACK(icnt).I = II;
PRTCL.TRACK(icnt).J = JJ;






return