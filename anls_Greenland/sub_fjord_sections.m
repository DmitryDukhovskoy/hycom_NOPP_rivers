function SGM = xsection_ts_fjords(HH,LON,LAT);
% Greenland fjords xsection
% for ARCc0.08
% topo 11

SGM(1).Name='WGr1';
SGM(1).IJ=[ 600        1010
         617         998
         625         993
         636         986];
SGM(2).Name='WGr2';
SGM(2).IJ=[583   769
   596   758
   607   745
   616   738];
SGM(3).Name='WGr3';
SGM(3).IJ=[560   705
   579   703
   602   694];
SGM(4).Name='EGr1';
SGM(4).IJ=[ 661   547
   668   541
   677   532];
SGM(5).Name='EGr2';
SGM(5).IJ=[   843   688
   844   683
   846   678
   846   670
   848   665
   851   661
   860   655
   871   654
   885   650];
SGM(6).Name='EGr3';
SGM(6).IJ=[852   721
   873   699
   881   694];
SGM(7).Name='EGr4';
SGM(7).IJ=[895   888
   908   892
   915   914
   923   917];

SGM(8).Name='EGr5';
SGM(8).IJ=[649 498;
	   661 488];
SGM(9).Name='WGr4';
SGM(9).IJ=[526 530;
	   543 530];

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