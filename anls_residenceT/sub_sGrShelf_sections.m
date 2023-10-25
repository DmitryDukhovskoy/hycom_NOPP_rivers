% South east and Southwest Greenland shelf sections
%
function SGM = sub_sGrShelf_sections(HH,LON,LAT);

SGM(1).Name='WGr1';
SGM(1).IJ=[510 556;
	   547 556];

SGM(2).Name='WGr2';
SGM(2).IJ=[520   459
   562   459];

SGM(3).Name='EGr1';
SGM(3).IJ=[630 473;
	   680 473];

SGM(4).Name='EGr2';
SGM(4).IJ=[617 438;
    660 438];



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
