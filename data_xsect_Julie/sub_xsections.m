% Section indices
% for T/S vertical cross-Arctic sections
function XSCT = sub_xsections(xname,res);
% these are indices for 0.08 - grid
switch(xname);
 case('BerNatl');
% 0.08 indices:
  IJs=[645, 1916;
       934, 1242;
       1038, 855;
       1020, 495;
       815, 179];
  ct1=-2;
  ct2=5;
  xl1=1;
  xl2=9450;
  yl1=-4500;
  yl2=0;

 case('BeaufNatl');
  IJs=[479, 1687;
       928, 1250
       1144, 272];
  ct1=-2;
  ct2=5;

 case('BeaufIcel');
  IJs=[479, 1687;
       928, 1250;
       1046 848;
       927 539];
  xl1=1;
  xl2=6735;
  yl1=-4500;
  yl2=0;
%  cs1=30;
%  cs2=35;  
  cs1=3.36;  % log scale S=28.79
  cs2=3.56;  % log scale S=35.16
  ct1=-2;
  ct2=5;

 case('LaptSea');
  IJs=[1096 1446
       1170 1546
       1170 1784];
  xl1=650;
  xl2=1420;
  yl1=-100;
  yl2=0;
  cs1=10;
  cs2=34;

 case('BafnFram');
  IJs=[  558         998
         497         644
         505         343
         638         140
         977         473
         960         622
        1041         882
         995        1044];
  xl1=1;
  xl2=10590;
  yl1=-5000;
  yl2=0;
  cs1=3.36;  % log scale S=28.79
  cs2=3.56;  % log scale S=35.16
  ct1=-2;
  ct2=5;

 case('BeringNorw');
%section specified by Wilbert - along 170W - 10 E
% Vertical transects of T and S
% Across the Arctic, Bering-to-North Pole-to-Nordic Seas along 
% the 170W/10E meridians, bounded by 65N (or the Norwegian coast).
  IJs= [657,    1945;
        679,    1849;
        722,    1717;
        809,    1489;
        901,    1305;
        1058,    981;
        1136,    771;
        1182,    632;
        1206,    547];

  xl1=1;
  xl2=7300;
  yl1=-5000;
  yl2=0;
  cs1=3.36;  % log scale S=28.79
  cs2=3.56;  % log scale S=35.16
  ct1=-2;
  ct2=5;

end

if res == 0.04
  IJs=IJs*2-1; % convert to 0.04 grid
end


% Find grid segments to connect large segments:
IIa=IJs(:,1);
JJa=IJs(:,2);

nij=size(IJs,1);
IIs=[];
JJs=[];
for ii=1:nij-1
  i1=IJs(ii,1);
  i2=IJs(ii+1,1);
  j1=IJs(ii,2);
  j2=IJs(ii+1,2);
  [I,J]=sub_xsct_indx(i1,j1,i2,j2);
  if size(I,1)==1;
    I=I';
    J=J';
  end

  IIs=[IIs;I];
  JJs=[JJs;J];
end;

IJs=[IIs,JJs];

XSCT.IJs = IJs;
XSCT.xl1 = xl1;
XSCT.xl2 = xl2;
XSCT.yl1 = yl1;
XSCT.yl2 = yl2;
XSCT.cs1 = cs1;
XSCT.cs2 = cs2;
XSCT.ct1 = ct1;
XSCT.ct2 = ct2;
XSCT.Isgm = IIa;
XSCT.Jsgm = JJa;

return 



