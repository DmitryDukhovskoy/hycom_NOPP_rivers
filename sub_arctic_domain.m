function ARC = sub_arctic_domain(HH,hmin);
% Prepare indices for Arctic domain
% hmin <0 
% exclude shallow regions with Hbtm>hmin
IJ = [   675        1920
         848        1899
        1067        1867
        1249        1769
        1383        1582
        1594        1241
        1579         647
        1501         587
        1370         664
        1276         720
        1129         910
        1104         956
         926         963
         855         970
         755        1041
         737        1114
         583        1234
         523        1364
         407        1449
         301        1586
         479        1846
         623        1920];

[mm,nn]=size(HH);

[II,JJ]=meshgrid([1:nn],[1:mm]);
IN = inpolygon(II,JJ,IJ(:,1),IJ(:,2));
HH(~IN) = 0;
Iarc = find(HH<=hmin);
Inan = find(HH>hmin | isnan(HH));

ARC.II   = IJ(:,1);
ARC.JJ   = IJ(:,2);
ARC.IN   = Iarc;
ARC.Inan = Inan;

return