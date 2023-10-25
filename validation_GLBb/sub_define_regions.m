function AA = sub_define_regions(HH);
% Define "boxes" - regions for 
%  heat and FWcontent analysis
%  
[mm,nn]=size(HH);
[X,Y]=meshgrid((1:nn),(1:mm));

% Labrador Sea
IJs=[119, 331;
     221, 331;
     221, 197;
     119, 197];
IJs(end+1,:)=IJs(1,:);
IMSK = inpolygon(X,Y,IJs(:,1),IJs(:,2));
IN = find(IMSK==1);
AA(1).Name = 'LabrSea';
AA(1).IJs  = IJs;
AA(1).IN   = IN;

% Subpolar Gyre -
% ~ in the center of the cold blob
IJs=[300, 30;
     300, 208;
     418, 208;
     418, 30];
IJs(end+1,:)=IJs(1,:);
IMSK = inpolygon(X,Y,IJs(:,1),IJs(:,2));
IN = find(IMSK==1);
AA(2).Name = 'SubpolarG';
AA(2).IJs  = IJs;
AA(2).IN   = IN;

% Greenland Gyre
IJs=[670, 560;
     670, 672;
     760, 672;
     760, 560];
IJs(end+1,:)=IJs(1,:);
IMSK = inpolygon(X,Y,IJs(:,1),IJs(:,2));
IN = find(IMSK==1);
AA(3).Name = 'GreenlandG';
AA(3).IJs  = IJs;
AA(3).IN   = IN;

% Irminger Sea
IJs=[330, 307;
     440, 307;
     440, 211;
     330, 211];
IJs(end+1,:)=IJs(1,:);
IMSK = inpolygon(X,Y,IJs(:,1),IJs(:,2));
IN = find(IMSK==1);
AA(4).Name = 'IrmingerSea';
AA(4).IJs  = IJs;
AA(4).IN   = IN;



return