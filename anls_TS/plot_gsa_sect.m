% Plot GSA sections 
% bbbbased on Belkin et al., 1998 
% also add FW pathway
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

s_fig=1;
regn = 'ARCc0.08';
expt = 110; % no Greenland runoff  
%expt = 112;  % Greenland runoff

fprintf('expt %3.3i\n\n',expt);


pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);

%ftopo = sprintf('%s/depth_%s_09.nc',pthtopo,regn); % old ARCc0.08
ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
Lmsk=HH*0;
Lmsk(HH<0)=1;

pmsk=[0,0,0; 1,1,1];

%[DX,DY]=sub_dx_dy(LON,LAT);
%Acell=DX.*DY; % Grid cell area, m2
clear SCT
SCT = sub_gsa_sections;
nst = length(SCT);

for ik=1:nst
  xy=SCT(ik).XY;
  IJ=[];
  for k=1:2
    x0=xy(k,1);
    y0=xy(k,2);
    D=distance_spheric_coord(LAT,LON,y0,x0);
    [j,i]=find(D==min(min(D)),1);
    IJ(k,1)=i;
    IJ(k,2)=j;
  end
  SCT(ik).IJ=IJ;
end

for ik=1:nst
  IJ=SCT(ik).IJ;
  i1=IJ(1,1);
  i2=IJ(2,1);
  j1=IJ(1,2);
  j2=IJ(2,2);
  [I,J]=sub_xsct_indx(i1,j1,i2,j2);
  SCT(ik).IIs=I;
  SCT(ik).JJs=J;
  SCT(ik).IND=sub2ind(size(HH),J,I);
end

figure(10); clf;
pcolor(Lmsk); shading flat;
colormap(pmsk); 
freezeColors;

hold on;
contour(HH,[-500:100:-10],'Color',[0.7 0.7 0.7]);
contour(HH,[-10000:1000:-900],'Color',[0.5 0.5 0.5]);
axis('equal');

for ik=1:nst
  II=SCT(ik).IIs;
  JJ=SCT(ik).JJs;
  plot(II,JJ,'.-','Color',[0 0.4 0.7]);
%  text(II(1),JJ(1),sprintf('%2.2i',ik),'Fontsize',12);
end

set(gca,'xlim',[350 1300],...
	'ylim',[10 1200],...
	'xtick',[],...
	'ytick',[]);

% Plot Greenland FW pathway
% see anls_UV/tracer_psthways_meanUV.m
npth=1;
pthm2 =sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',110);
fmpth =sprintf('%strcr_pathway_%2.2i.mat',pthm2,npth);
load(fmpth);
IIs=PTS.IJpathway(:,1);
JJs=PTS.IJpathway(:,2);
Dst=PTS.Dist;

plot(IIs,JJs,'r.','Markersize',12);
% Plot original nodes:
%plot(X,Y,'k^','Markersize',5,...
%          'MarkerFaceColor','k');
for ik=1:11
  xx=(ik-1)*1000;
  D=abs(Dst-xx);
  i0=find(D==min(D),1);
  plot(IIs(i0),JJs(i0),'b.','Markersize',26);
  stl=sprintf('%i',xx);
%  text(IIs(i0),JJs(i0),stl,'Fontsize',12);
end



btx='plot_gsa_sect.m';
bottom_text(btx,'pwd',1);

if s_fig==1
  fgnm=sprintf('%sgsa_sections.png',pthfig);
  fprintf('Saving %s\n',fgnm);
  print('-dpng','-r250',fgnm);
end


  



