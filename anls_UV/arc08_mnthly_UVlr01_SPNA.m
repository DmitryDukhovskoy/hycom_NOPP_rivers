% subsample monthly fileds 1st surface layer UV into SPNA
% for Igor Yashayaev analysis
% see mnthly_arc08_UV.m
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

YR1 = 2000;
YR2 = 2016;

s_mat = 1; % =0 - do not save mat file
           % =1 - save tracer and overwrite existing mat
	   % =2 - skip months where mat file exist

if s_mat==0,
  fprintf('Mat file is not created\n');
elseif s_mat == 1
  fprintf('Mat file will be saved, old mat file will be overridden\n');
elseif s_mat == 2
  fprintf('Extraction is skipped for months where old mat files exist\n');
end

% Surface layer 01
plr   = 1;  % layer to calculate mean U

regn = 'ARCc0.08';
expt = 110;
rg = 9806;

pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);


fprintf('Monthly mean UV %i-%i in the layer: %i\n',YR1,YR2,plr);

%figure(1); clf;
%set(gcf,'Visible','off');


% SPNA region:
IJ=[     365         766
         365          45
        1236          45
        1236         766];

ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);

IG=IJ(:,1);
JG=IJ(:,2);
[IDX,JDX]=meshgrid([1:nn],[1:mm]);
IN=inpolygon(IDX,JDX,IG,JG);
IRG=find(IN==1);


hmsk=HH;
hmsk(HH<0)=nan;

if ~exist('DX','var')
  [DX,DY]=sub_dx_dy(LON,LAT);
end

i1=min(IJ(:,1));
i2=max(IJ(:,1));
j1=min(IJ(:,2));
j2=max(IJ(:,2));
HHs=HH(j1:j2,i1:i2);

UV.Bath_m=HHs;
UV.LON=LON(j1:j2,i1:i2);
UV.LAT=LAT(j1:j2,i1:i2);

cc=0;
Umn=[];
Vmn=[];
for iyr=YR1:YR2
  fmatin = sprintf('%smnthUV_lr%2.2i_%i.mat',pthmat,plr,iyr);

		fprintf('Loading %s\n',fmatin);
		load(fmatin);

% Average
  for imo=1:12
    U=meanUV(imo).U;
    V=meanUV(imo).V;
 
    cc=cc+1;
    if isempty(Umn)
      Umn=U;
      Vmn=V;
    else
      Umn=Umn+U;
      Vmn=Vmn+V;
    end 
%    UV(imo).U=U(j1:j2,i1:i2);
%    UV(imo).V=V(j1:j2,i1:i2);
    fprintf('cc=%i\n',cc);
  end

%  fmatout = sprintf('%sarc08_mnthUVspna50m_%i.mat',pthmat,iyr);
%  fprintf('Saving %s\n',fmatout);
%  save(fmatout,'UV');
end
Umn=Umn/cc;
Vmn=Vmn/cc;
UV.Umean=Umn(j1:j2,i1:i2);
UV.Vmean=Vmn(j1:j2,i1:i2);

fprintf('Max/min v = %5.2f %5.2f, max/min v=  %5.2f/%5.2f\n',max(max(Umn)),...
         min(min(Umn)), max(max(Vmn)), min(min(Vmn)));

fmatout = sprintf('%sarc08_mnthUVspna_lr%2.2i_%i-%i.mat',pthmat,plr,YR1,YR2);
fprintf('Saving %s\n',fmatout);
save(fmatout,'UV');

 



