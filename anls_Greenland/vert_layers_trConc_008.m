% Compute mean concentration
% of the tracer by layers
% within the specified domain IN
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

nTr   = 1   % tracer to plot
plr   = 0;  % layer to integrate, =0 <- over whole depth, use this
z0    = -500; % min depth

s_mat = 2; % =1 - save mat file flag
           % =2 - load saved fields

TrMn=1e-12; % Threshold value to plot
rg = 9806;
fprintf('Tracer #: %i, Threhold value: %8.5d\n',nTr,TrMn);


rg = 9806;
YRPLT=[];
cc=0;
for yr=1993:2016
  for im=1:12
    iday=datenum(yr,im,15)-datenum(yr,1,1)+1;
    cc=cc+1;
    YRPLT(cc,1)=yr;
    YRPLT(cc,2)=iday;
  end
end



nyy = size(YRPLT,1);

% Experiments:
pthfig = sprintf('/Net/mars/ddmitry/hycom/%s/%3.3i/fig_trac/',regn,expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc0.08/%3.3i/data_mat/',expt);

fmat = sprintf('%s%s_%3.3i_LrAvrg_Tr%2.2i_z%im_%i_%i.mat',...
		   pthmat,regn,expt,nTr,abs(z0),YRPLT(1,1),YRPLT(end,1));


%ftopo = sprintf('%s/depth_%s_17DD.nc',pthtopo,regn); % 
ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[DX,DY]=sub_dx_dy(LON,LAT);
Acell=DX.*DY; % Grid cell area, m2
MVOL = [];
HHm=HH;
HHm(HHm>z0)=100;

if s_mat<2
  for ip=1:nyy;
    yr   = YRPLT(ip,1);
    iday = YRPLT(ip,2);
    dnmb = datenum(yr,1,1)+iday-1;
    DV   = datevec(dnmb);

    %pthbin = sprintf('/nexsan/hycom/ARCc0.08_011/data/%i/',yr);  
    pthbin = sprintf('/nexsan/archive/%s_%3.3i/data/%i/',regn,expt,yr);
    fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
    finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);

    fprintf('\n:::: Analyzing %i/%2.2i/%2.2i   ::::\n\n',DV(1:3));

    tic;
    fprintf('Computing layer-averaged concentration Tracer %i ...\n',nTr);

    if plr>0
      [F,n,m,l] = read_hycom(fina,finb,'tracer','r_tracer',nTr,'r_layer',plr);
      F(F>1e6)=nan;
      Tr=F(plr,:,:);

      [F,n,m,l] = read_hycom(fina,finb,'thknss','r_layer',plr);
    %  F=squeeze(F(plr,:,:));
      F=F./rg;
      F(F>1e10)=0;
      F(F<1e-2)=0;
      dH=F; 

    else
      [F,n,m,l] = read_hycom(fina,finb,'tracer','r_tracer',nTr);
      F(F>1e6)=nan;
      Tr=F;

      [F,n,m,l] = read_hycom(fina,finb,'thknss');
      F=F./rg;
      F(F>1e10)=0;
      F(F<1e-2)=0;
      dH=F; 

      dmm = sub_layer_avrg_trConc(dH,Tr,HHm,Acell);
      
      TRI.Tracer                  = nTr;
      TRI.avrg_TrConc_kg_m3(:,ip) = dmm.Vert_ConcMean_kg_m3;
      TRI.Mean_LThkn_m(:,ip)      = dmm.Mean_LThkn_m;
      TRI.TM(ip)                  = dnmb;
      
%keyboard

      if s_mat==1 & (mod(ip,6)==0 | ip==nyy)
%	fmat = sprintf('%s%s_%3.3i_VolIntgr_Tr%2.2i_%i%2.2i%2.2i.mat',...
%		       pthmat,regn,expt,nTr,DV(1:3));
	fprintf('Saving %s\n',fmat);
	save(fmat,'TRI');
      end

      fprintf('1 record processed %8.6fm\n\n',toc/60);
    end % time

  end; % time loop

else
  fprintf('Loading %s\n',fmat);
  load(fmat);
    
end

TM = TRI.TM;
dzm= TRI.Mean_LThkn_m;
Ctr= TRI.avrg_TrConc_kg_m3;
Ctr(Ctr<TrMn)=nan;
lgC= log(Ctr);

ZM = -cumsum(dzm);
yl1 = round(min(min(ZM)));

DV = datevec(TM);
nrc = length(TM);
yr1= DV(1,1);
yr2= DV(end,1);
clear yrs
for ik=1:nrc
  yrs(ik,1)=DV(ik,1)+(DV(ik,2)-1)/12;
end

[YRS,dmm]=meshgrid(yrs,ZM(:,1));



nint = 320;
c1 = -4;
c2 = -1;
%CMP = create_colormap2_3(nint,c1,c2);
CMP = colormap_sclr2(nint,c1,c2);
cmp = CMP.colormap;
for ik=1:20
  cmp(ik,:) = [1 1 1];
end
cmp = smooth_colormap(cmp,20);
cmp = smooth_colormap(cmp,20);
cmp(1,:) = [1 1 1];
cnt = CMP.intervals;


figure(1); clf;
axes('Position',[0.1 0.5 0.84 0.43]);
pcolor(YRS,ZM,lgC); shading interp;
colormap(cmp);
caxis([c1 c2]);
set(gca,'tickdir','out',...
	'xtick',[1990:2016],...
	'xlim',[yr1 yr2],...
	'ylim',[-4000 0],...
	'ytick',[-4500:500:0]);


nm = sprintf('%s-%3.3i',regn,expt);
stl = sprintf('%s, Conc. mean layers, log(Tr) %i, kg/m3',nm,nTr);
title(stl);

hc=colorbar('Location','SouthOutside');
set(hc,'position',[0.1 0.41 0.84 0.012],...
    'Ticklength',0.01,'Fontsize',12);

btx = 'vert_layers_trConc_008.m';
bottom_text(btx,'pwd',1,'position',[0.08 0.33 0.8 0.1]);




