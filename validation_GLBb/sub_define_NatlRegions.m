function AA = sub_define_NatlRegions(HH,LON,LAT);
% Define "boxes" - regions for 
%  heat and FWcontent analysis
%  Use same boxes as for EN4 analysis
%addpath /home/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers/EN4_anls
fprintf('Defining regions, same as EN4 ...\n');
[mm,nn]=size(HH);
[X,Y]=meshgrid((1:nn),(1:mm));

BB = sub_regions_EN4;
nb = length(BB);

for ii = 1:nb
  nm = BB(ii).Name;
  XY = BB(ii).XY;

  clear IJs
  for jj=1:4
    x1=XY(jj,1);
    y1=XY(jj,2);
    
    dd = distance_spheric_coord(LAT,LON,y1,x1);
    [j0,i0] = find(dd==min(min(dd)));
    if dd>15000,
      fprintf('Searching for X=%6.3fE Y=%6.3fN\n',x1,y1);
      error('EN4 box is outside HYCOM NAtl subdomain');
    end
    
    
    IJs(jj,1)=i0;
    IJs(jj,2)=j0;
  end    
  
  IJs(end+1,:)=IJs(1,:);
  
  IMSK = inpolygon(X,Y,IJs(:,1),IJs(:,2));
  IN = find(IMSK==1);
  
  AA(ii).Name = nm;
  AA(ii).IJs  = IJs;
  AA(ii).IN   = IN;

end



return