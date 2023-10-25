function RG = regions;
% Define regions: Canada Basin, Eurasian Basin

% Canada Basin
CB = [  1022        1526
         793        1669
         653        1672
         502        1666
         465        1615
         467        1500
         615        1327
         745        1206
         851        1185];

CB(end+1,:) = CB(1,:);

% Eurasian Basin:
EB = [        1055        1526
         862        1116
         963         986
        1070         998
        1156        1244
        1194        1491
        1149        1568];
EB(end,:) = EB(1,:);

RG(1).Title = 'Canada Basin';
RG(1).IJ    = CB;

RG(2).Title = 'Eurasian Basin';
RG(2).IJ    = EB;


return