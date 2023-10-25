% Get indices from CFSv2 grid
% to interpolate onto HYCOM ARCc0.72 grid
% Use closest neighbour approach
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08/;
addpath /usr/people/ddmitry/codes/MyMatlab/;
startup;

clear
close all

year=2011;
mnth=12;
dd=1;
%pthdat=(['/Net/Movies0/ddmitry/vector_winds_ccmp_level30/',int2str(year),'/']);
pthd = '/Net/kronos/ddmitry/ncep_cfsv2/';
pth72='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.72/topo_grid/';
%pth8='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat='/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.72/data_mat/';

%ffout=[pth8,'cfsv2_gridindx_arc008all.mat'];
%ffoutA=[pth8,'cfsv2_gridindx_arc008TMP.mat'];
fmat=sprintf('%scfs_gridindx_arc072_nghb.mat',pthmat);


date_s=sprintf('%4i%2.2i%2.2i',year,mnth,dd)
fnm = sprintf('cfsv2-sec2_%4.4i_mon_uv-10m.nc',year);
fp = [pthd,fnm];

s_mat = 1;	% =1 - save final mat file; =2 save temporary and
                % final; =0 - do not save mat files   

if s_mat==0
  fprintf('output mat files will not be created\n\n');
end


X = nc_varget(fp,'Longitude');
X(X>180)=X(X>180)-360;
Y = nc_varget(fp,'Latitude');
nc=length(X);
mc=length(Y);



% Topo for ARCc0.72
fsv=[pth72,'new_bath072.mat'];
load(fsv);
LN=elon;
LT=alat;
HH=hnew;
HH(isnan(HH))=10;
clear hnew elon alat;
[mm,nn]=size(HH);

ni = mm*nn;
INDX.dim_rows = mc;
INDX.dim_colms= nc;

for ik=1:ni
  if mod(ik,1000)==0
    fprintf('Searching for indices, %5.1f done ...\n',(ik/ni)*100);
  end
  
  [j0,i0]=ind2sub(size(HH),ik);
  x0=LN(ik);
  y0=LT(ik);
  dy = abs(Y-y0);
  dx = abs(X-x0);
  jh0 = find(dy==min(dy),1);
  ih0 = find(dx==min(dx),1);
  IH0 = sub2ind([mc,nc],jh0,ih0);
  INDX.II(j0,i0) = IH0;
end

if s_mat==1
  fprintf('Saving output %s\n',fmat);
  save(fmat,'INDX');
end


