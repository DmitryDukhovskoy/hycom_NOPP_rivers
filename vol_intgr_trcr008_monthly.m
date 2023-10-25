% Compute volume integrated mass of the tracer
% "monthly" - 15th day of every month, not monthly mean
% within the specified domain IN
% Use only deep basin - deeper than 500m
% If IN is empty - the whole domain
% nTr - tracer # that is integrated
% Acell - grid cell area, m2
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

regn = 'ARCc0.08';
expt = 110;  

nTr   = 1;   % tracer to plot
s_mat = 1; % =1 - save mat file flag
           % =2 - load saved fields
YR1 = 2015;
YR2 = 2015;
h0  = -500; % cutoff depth, m

rg = 9806;


% Experiments:
pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc0.08/%3.3i/data_mat/',expt);


%ftopo = sprintf('%s/depth_%s_17DD.nc',pthtopo,regn); % 
ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2
MVOL = [];
Hmsk = HH;
Hmsk(HH>h0)=100;

fprintf('Tracer: %i, Depths<%4.1fm, Layer integrated, %i-%i\n',nTr,h0,YR1,YR2);
if s_mat==1
  fprintf('Mat file will be saved\n');
else
  fprintf('Mat file WILL NOT be saved\n');
end  

for iyr=YR1:YR2
  yr   = iyr;
  fmat = sprintf('%s%3.3i_VrtLrMass_mo_Tr%2.2i_%i.mat',...
		 pthmat,expt,nTr,iyr);
  cc=0;
  TRI.Title = sprintf('%s Tracer Mass integrated over domain by layers',regn);
  for imo=1:12
    dnmb = datenum(yr,imo,15);
    DV   = datevec(dnmb);
    iday = dnmb-datenum(yr,1,1);
  
%pthbin = sprintf('/nexsan/hycom/ARCc0.08_011/data/%i/',yr);  
    pthbin = sprintf('/nexsan/archive/%s_%3.3i/data/%i/',regn,expt,yr);
    fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
    finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);

    fprintf('\n:::: Analyzing %i/%2.2i/%2.2i   ::::\n\n',DV(1:3));
    cc=cc+1;
    
    tic;
  
%    fprintf('Computing vol-integrated Tracer mass ...\n');
    [F,n,m,l] = read_hycom(fina,finb,'tracer','r_tracer',nTr);
    F(F>1e6)=nan;
    Tr=F;

    [F,n,m,l] = read_hycom(fina,finb,'thknss');
    F=F./rg;
    F(F>1e10)=0;
    F(F<1e-2)=0;
    dH=F; 

    tri = sub_intgr_tr(dH,Tr,Hmsk,Acell);
    TRI.TM(cc) = dnmb;
    TRI.Vert_LrMass_kg(:,cc) = tri.Vert_LrMass_kg;
    TRI.Mean_LThkn(:,cc) = tri.Mean_LThkn;
    
    fprintf('1 day processed %6.2fmin\n\n',toc/60);
    
  end
  
  if s_mat==1 
    fprintf('Saving %s\n',fmat);
    save(fmat,'TRI');
  end

end; % time loop

