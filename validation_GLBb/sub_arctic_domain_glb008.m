function ARC = sub_arctic_domain(HH,hmin);
% Prepare indices for Arctic domain
% hmin <0 
% exclude shallow regions with Hbtm>hmin
IJ = [   458        1771
         602        1771
        1031        1630
        1308        1271
        1363        1131
        1365         513
        1271         437
        1057         572
         918         760
         905         800
         868         802
         721         822
         678         797
         549         903
         514         971
         362        1075
         312        1209
         253        1263
         150        1395
         159        1521
         314        1720
         416        1771];

[mm,nn]=size(HH);

[II,JJ]=meshgrid([1:nn],[1:mm]);
IN = inpolygon(II,JJ,IJ(:,1),IJ(:,2));
HH(~IN) = 100;
Iarc = find(HH<=hmin);
Inan = find(HH>hmin | isnan(HH));
%keyboard

ARC.II   = IJ(:,1);
ARC.JJ   = IJ(:,2);
ARC.IN   = Iarc;
ARC.Inan = Inan;

return