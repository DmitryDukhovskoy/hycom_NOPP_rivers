function SGM = xsection_ts_fjords(HH,LON,LAT);
% Greenland fjords xsection
% for ARCc0.08
% topo 11

SGM(1).Name='WGr1';
SGM(1).IJ=[1202        2019
        1235        1994
        1268        1974]; 
    
SGM(2).Name='WGr2';
SGM(2).IJ=[1169        1538
        1196        1506
        1219        1483
	1218        1480
        1218        1469
        1220        1463
        1232        1449];
SGM(3).Name='WGr3';
SGM(3).IJ=[1124        1415
        1156        1402
        1199        1385
        1218        1378];
SGM(4).Name='EGr1';
SGM(4).IJ=[1323        1092
        1334        1080
        1345        1067];
SGM(5).Name='EGr2';
SGM(5).IJ=[1664        1384
        1678        1382
	1684        1375   
        1688        1363
        1694        1327
        1721        1309
        1761        1300];
SGM(6).Name='EGr3';
SGM(6).IJ=[1704        1442
        1719        1424
        1760        1381];
SGM(7).Name='EGr4';
SGM(7).IJ=[1791        1780
        1815        1784
        1821        1796
        1825        1810
        1827        1824
        1840        1829];

nsgm = length(SGM);

for k = 1:nsgm
  IJs = SGM(k).IJ;
  nij = size(IJs,1);
  xnm = SGM(k).Name;
  
  if strncmp(xnm,'WGr',3),
    xdr=-1;
  else
    xdr=1;
  end

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

    if ~isempty(IIs) & I(1)==IIs(end),
      I=I(2:end);
      J=J(2:end);
    end
    
    IIs=[IIs;I];
    JJs=[JJs;J];
  end;

  IJs=[IIs,JJs];

  nS=length(IIs);
  clear Xl Yl
  if xdr>0
    i1=1;
    i2=nS;
    di=1;
  else
    i1=nS;
    i2=1;
    di=-1;
  end
%  keyboard

  cii=0;
  clear dx
  for ii=i1:di:i2 % distance from the coast, land=0m
    cii=cii+1;
    i0=IJs(ii,1);
    j0=IJs(ii,2);
    Xl(ii,1)=LON(j0,i0);
    Yl(ii,1)=LAT(j0,i0);
    if ii==i1
      dx(ii,1)=0;
    else
      dx(ii,1)=distance_spheric_coord(Yl(ii-di),Xl(ii-di),Yl(ii),Xl(ii));
    end
  end;
  if xdr>0
    dst=cumsum(dx);
  else
    dst = cumsum(dx,'reverse');
  end
  
  INDs=sub2ind(size(HH),JJs,IIs);
  
%  SGM(k).IJs   = IJs;
  SGM(k).IIs   = IIs;
  SGM(k).JJs   = JJs;
  SGM(k).Xcrd  = Xl;
  SGM(k).Ycrd  = Yl;
  SGM(k).INDs  = INDs;
  SGM(k).dx_m  = dx;
  SGM(k).dist_m= dst;
  
end
  
SGM = rmfield(SGM,'IJ');


return