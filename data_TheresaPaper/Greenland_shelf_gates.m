% Coordinates of the sections (gates)
% on Gr shelf for calculating fluxes
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

YR1=1993;
YR2=2016;
zz1=0;
zz2=-50;

av=0; % =1 - average over specofoed years, =0 - plot individual years
plr =1;  % U from plr layer
rg = 9806;
rhow=1027; 
hgg=1e20;

regn = 'ARCc0.08';
expt = 110;
pthfig  = sprintf('/nexsan/people/ddmitry/Net_ocean/hycom/ARCc0.08/%3.3i/fig_meanUV/',expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat =sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
pthmat2 =sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat2/',expt);
pthout='/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_theresa/';

ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);

GC = sub_greenl_isobath(HH,LON,LAT);
Ig=GC.cntr_Iindx;
Jg=GC.cntr_Jindx;


% Specify sections:
% Greenlaaand    shelf OSNAP mooring line
SCT(1).Name='DenmarkStr';
SCT(1).IJ=[788 610; 847 536];
SCT(2).Name='EastGrShelf_OSNAP';
SCT(2).IJ=[   613   425; 671   405];
SCT(3).Name='DavisStr';
SCT(3).IJ=[471 673; 556 652];
SCT(4).Name='NaresStr';
SCT(4).IJ=[570 1037; 605 999];
SCT(5).Name='FramStr';
SCT(5).IJ=[932 960; 1072     960];
nsct=length(SCT);

for ip=1:nsct
  fprintf('Section: %s\n',SCT(ip).Name);
  IJs=SCT(ip).IJ;
  [IIs,JJs]=sub_xsct_indx(IJs(1,1),IJs(1,2),IJs(2,1),IJs(2,2));
  SCT(ip).I=IIs;
  SCT(ip).J=JJs;
end

% Find normal at every small segments (+ northward, + eastward)
% Find intersection with GR contour
for ip=1:nsct
  IIs=SCT(ip).I;
  JJs=SCT(ip).J;
  ni=length(IIs);
  ist=IIs(1);
  jst=JJs(1);
  
  clear Nrm Hbtm dL Dst Xlon Ylat
  for isg=1:ni-1
    i0 = IIs(isg);
    i1 = IIs(isg+1);
    j0 = JJs(isg);
    j1 = JJs(isg+1);

    if i1==i0 % Y section, X norm ->
      Nrm(isg,1)=1;
      Nrm(isg,2)=0;  
    else    % Y norm ^
      Nrm(isg,1)=0;
      Nrm(isg,2)=1;
    end
    Hbtm(isg,1)=HH(j0,i0);
% Segments length:    
    dl=distance_spheric_coord(LAT(j0,i0),LON(j0,i0),LAT(j1,i1),LON(j1,i1));
    dL(isg,1)=dl;
    dst=distance_spheric_coord(LAT(j0,i0),LON(j0,i0),...
			       LAT(jst,ist),LON(jst,ist));
    Dst(isg,1)=dst;
    Xlon(isg,1)=LON(j0,i0);
    Ylat(isg,1)=LAT(j0,i0);
  end
  Hbtm(ni)=HH(j1,i1);
  Dst(ni)=distance_spheric_coord(LAT(j1,i1),LON(j1,i1),...
			       LAT(jst,ist),LON(jst,ist));
  dL(ni)=dl;
  
  SCT(ip).Nrm=Nrm;
  SCT(ip).Hbottom=Hbtm;
% Intersection pnt:
  dd=sqrt((IIs-Ig).^2+(JJs-Jg).^2);
  [jm,im]=find(dd==min(min(dd)),1);
  imn=IIs(im);
  jmn=JJs(im);
  
  SCT(ip).GrCntr_intrcp=im;
  SCT(ip).GrCntr_I=imn;
  SCT(ip).GrCntr_J=jmn;
  SCT(ip).Segm_dL=dL;
  SCT(ip).Dist_origin=Dst;
  SCT(ip).longt=Xlon;
  SCT(ip).latit=Ylat;
  
end


f_map=0;
if f_map==1
  fn=1;
  sub_plot_Greenl_contour(HH,LON,LAT,fn,GC);

  for ip=1:nsct
%    IJ=SCT(ip).IJ;
%    plot([IJ(1,1) IJ(2,1)],[IJ(1,2) IJ(2,2)],...
%	 'Linewidth',2.5,'Color',[1. 0.6 0]);
    IIs=SCT(ip).I;
    JJs=SCT(ip).J;
    plot(IIs,JJs,'-',...
	 'Linewidth',2.5,'Color',[1. 0.6 0]);
  end
  bottom_text(btx,'pwd',1);
end

fmat=sprintf('%shycom008_GreenlandGates',pthout);
fprintf('Saving %s\n',fmat);
save(fmat,'SCT');
