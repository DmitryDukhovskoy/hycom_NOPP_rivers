function AA = sub_define_regionsAO(HH);
% Define "boxes" - Arctic Ocaen regions for 
%  heat and FWcontent analysis
%  
[mm,nn]=size(HH);
[X,Y]=meshgrid((1:nn),(1:mm));

% CAA
IJs=[50         790
     33        1333
     231       1341
     660        879
     438        790];
IJs(end+1,:)=IJs(1,:);
IMSK = inpolygon(X,Y,IJs(:,1),IJs(:,2));
IN = find(IMSK==1);
AA(1).Name = 'CAA';
AA(1).IJs  = IJs;
AA(1).IN   = IN;

% Canada Basin
IJs=[192        1383
     382        1563
     819        1400
     614         976];
IJs(end+1,:)=IJs(1,:);
IMSK = inpolygon(X,Y,IJs(:,1),IJs(:,2));
IN = find(IMSK==1);
AA(2).Name = 'Canada Basin';
AA(2).IJs  = IJs;
AA(2).IN   = IN;

% Eurasian Basin
IJs=[858       1412
    952        1422
    991        1319
    842         832
    620         965];
IJs(end+1,:)=IJs(1,:);
IMSK = inpolygon(X,Y,IJs(:,1),IJs(:,2));
IN = find(IMSK==1);
AA(3).Name = 'Eurasian Basin';
AA(3).IJs  = IJs;
AA(3).IN   = IN;

% East Euras. Shelf
IJs=[354        1731
     371        1551
     589        1541
     845        1435
     968        1430
     989        1299
    1175        1311
    1008        1694
     737        1756
     570        1738];
IJs(end+1,:)=IJs(1,:);
IMSK = inpolygon(X,Y,IJs(:,1),IJs(:,2));
IN = find(IMSK==1);
AA(4).Name = 'EastEuras Shelf';
AA(4).IJs  = IJs;
AA(4).IN   = IN;

% West Euras. Shelf
IJs=[1071        1323
     962        1190
     965        1025
     906         805
    1074         585
    1274         479
    1362         980
    1353        1148];
IJs(end+1,:)=IJs(1,:);
IMSK = inpolygon(X,Y,IJs(:,1),IJs(:,2));
IN = find(IMSK==1);
AA(5).Name = 'WestEuras Shelf';
AA(5).IJs  = IJs;
AA(5).IN   = IN;



return