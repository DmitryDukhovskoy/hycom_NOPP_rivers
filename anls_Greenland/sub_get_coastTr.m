   function TRC = sub_get_coastTr(HH,dH,Tr,mTr,LON,LAT,SGM,plr,dst,dd,dnmb);
% SGM - segment on land, 
% tracer is saved along the transect perpendicular
% to this segment
%
np = length(SGM.I);
cc = 0;
ig=5;
Io = find(HH<0 & HH>-10);
[jxo,ixo]=ind2sub(size(HH),Io);
Il = find(HH>=0);
[jxl,ixl]=ind2sub(size(HH),Il);
Iout = SGM.Iout; 
Jout = SGM.Jout;

fprintf('Getting data along transects, np=%i ...\n',np);

for k=1:dd:np
%  fprintf('k=%i\n',k);
  cc=cc+1;
  i1=SGM.I(k);
  j1=SGM.J(k);
%  i2=SGM.I(k+ig);
%  j2=SGM.J(k+ig);
%  I1 = sub2ind(size(HH),j1,i1);
  
% Start from the coast line: find closest pnt
% First move to the ocean, than on land
  dmm = abs(sqrt((ixo-i1).^2+(jxo-j1).^2));
  Imm = find(dmm==min(dmm),1);
  jmm = jxo(Imm);
  imm = ixo(Imm);
  
  dmm = abs(sqrt((ixl-imm).^2+(jxl-jmm).^2));
  Ill = find(dmm==min(dmm),1);
  ill = ixl(Ill);
  jll = jxl(Ill);
  i1 = ill;
  j1 = jll;

% Find normal to the coast line
% directed offshore
%
% Use predefined rotation vector rtr
%  a=i2-i1;
%  b=j2-j1;
%  c=rtr*(a+i*b);
%  ec=c/sqrt(c*conj(c));
%  iE = round(i1+(dst+1)*real(ec));
%  jE = round(j1+(dst+1)*imag(ec));
%
% Use outside contour
  dmm = sqrt((Iout-i1).^2+(Jout-j1).^2);
  ib = find(dmm==min(dmm),1);
  i2 = Iout(ib);
  j2 = Jout(ib);
  a=i2-i1;
  b=j2-j1;
  c=a+i*b;
  ec=c/sqrt(c*conj(c));
  iE = round(i1+(dst+1)*real(ec));
  jE = round(j1+(dst+1)*imag(ec));

  [I,J]= sub_xsct_indx(i1,j1,iE,jE);
  Isct = sub2ind(size(HH),J,I);

  DL = distance_spheric_coord(LAT(Isct),LON(Isct),LAT(j1,i1),LON(j1,i1));
%keyboard  
%  if (cc==10), keyboard; end;

  Isct = Isct(1:dst);
  TRC.dnmb              = dnmb;
  TRC.Nmb_vLaeyr_avrgd  = plr;
  TRC.Dist_m(cc,1:dst)  = DL(1:dst);
  TRC.Dpth(cc,1:dst)    = HH(Isct);
  TRC.Indx(cc,1:dst)    = Isct;
  TRC.X(cc,1:dst)       = LON(Isct);
  TRC.Y(cc,1:dst)       = LAT(Isct);
  TRC.Tr(cc,1:dst)      = Tr(Isct);
  TRC.mass_Tr(cc,1:dst) = mTr(Isct);
end;  
  
   
return