% Plot vertical disctribution
% of Greenland tracer
% on South Gr shelves
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/seawater;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

nTr   = 1;
s_mat = 1; % =0 - extract data, no save; =1 extract & save, =2 - load saved
s_map = 0; % plot map with segments
%f_tracloc = 0; % =1 - add locations of tracer release to transect maps
rg=9806;  % convert pressure to depth, m
%fld0='salin'; %'temp' - background field

plr=0; % highlight this interface
btx = 'xsection_tracer_fjords008.m';


regn = 'ARCc0.08';
%expt = 110; % experiment without runoff
expt = 112;  % epxeriment with Greenland runoff and monthly Arctic rivers
pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
pthtopo= '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';


ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);


YRPLT=[];
cc=0;
YR1=2016;
YR2=2016;
for iyr=YR1:YR2
  for idd=1:7:365
%  for idd=45:45
    cc=cc+1;
    YRPLT(cc,1)=iyr;
    YRPLT(cc,2)=idd;
  end
end

np   = size(YRPLT,1);
SGM  = sub_sGrShelf_sections(HH,LON,LAT);
nsct = length(SGM);

for isct=1:nsct
		SCT(isct).IJ_indx = [SGM(isct).IIs, SGM(isct).JJs];
		SCT(isct).segm_dL = SGM(isct).dx_m;
		SCT(isct).dist_m = SGM(isct).dist_m;

  II=SGM(isct).INDs;
  Hb=HH(II);
  SCT(isct).Hbottom=Hb;
  SCT(isct).Name = SGM(isct).Name;
end


f_map_check=0; % quickly see the segments
if f_map_check>0
  figure(10); clf;
%  fn = 10;
%  sub_plot_bath(HH,LON,LAT,fn,domname);
  hold on
  contour(HH,[0 0],'k');
  contour(HH,[-1000:100:-100],'Color',[0.8 0.8 0.8]);
  contour(HH,[-5000:1000:-1000],'Color',[0.9 0.9 0.9]);
  for kk=1:nsct
    IIs = SGM(kk).IIs;
    JJs = SGM(kk).JJs;
    plot(IIs,JJs,'b.-');
    stx=sprintf('Sect %2.2i %s',kk,SGM(kk).Name);
    text(min(IIs),min(JJs),stx);
  end
  axis('equal');
  set(gca,'xlim',[500 960],...
	  'ylim',[400 1080]);
  bottom_text(btx,'pwd',1);
end

%fprintf('Section: %s, Saving fig: %i\n',xname,s_fig);
%keyboard
% Interpolate on Z
ZZf = [(0:-1:-10)';(-12:-2:-200)';(-205:-5:-500)';...
       (-510:-10:-2000)';(-2025:-25:-2500)';(-2600:-100:-5000)'];
kzz = length(ZZf);

dZf=abs(diff(ZZf));
ZMf = [];
for ik=1:kzz-1
  ZMf(ik,1)=ZZf(ik)+0.5*(ZZf(ik+1)-ZZf(ik));
end

for isct=1:nsct
  SCT(isct).ZZ=ZZf;
end

fmatout = sprintf('%shycom008_%3.3i_Trcr%2.2i_sGrShxsct_%4.4i-%4.4i.mat',...
                    pthmat,expt,nTr,YR1,YR2);

% Plot fields:
cc=0;  % 
ip1=1;
Nlast=0;
for ip=ip1:np
  yr=YRPLT(ip,1);
  iday=YRPLT(ip,2);
  pthbin = sprintf('/nexsan/archive/%s_%3.3i/data/%i/',regn,expt,yr)
  
  if expt==112
    pthbin = sprintf('/nexsan/hycom/ARCc0.08_112/data/%i/',yr);
  end
    
  fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
  finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);

  
  if ~exist(fina,'file');
    fprintf('Not found: %s\n\n',fina);
    continue;
  end
  
  cc=cc+1;
  dnmb=datenum(yr,1,1)+iday-1;
  DV=datevec(dnmb);

  fprintf('%4.4i_%2.2i_%2.2i: %s\n',DV(1:3),fina);

  tic;
%  [F,n,nlev] = read_hycom(fina,finb,'temp');
  [F,n,m,l] = read_hycom(fina,finb,'tracer','r_tracer',nTr);
  F(F>1e6)=nan;
  F(F<0)=0;
  TR=F;  % kg/m3


  fld='thknss';
  [F,n,m,l] = read_hycom(fina,finb,fld);
  F(F>1e18)=nan;
  F=F/rg;
  F(F<1e-2)=0;
  dH = F;

% All sections are along X - need only V components - normal to xsct
% These are total velocities    
% barotropic not needed
		[F,n,m,nlr] = read_hycom(fina,finb,'v-vel.');
  hgg=1e20;
		F(F>hgg)=nan;
		V=F;

%		if arch_mean==0
%				fprintf('STOP: Add U,V barotropic to read from archv !!! \n');
%				error('Add u, v barotropic to the code');
%		end

		[ZMh,ZZh] = sub_zz_zm(fina,finb,HH);
		dH=abs(diff(ZZh,1));

  for isct=1:nsct
    INDs = SGM(isct).INDs;
    xnm  = SGM(isct).Name;
    Tr   = squeeze(TR(:,INDs));
    Un   = squeeze(V(:,INDs));
    dHn  = squeeze(dH(:,INDs));

    Is = SGM(isct).IIs;
    Js = SGM(isct).JJs;
    nns = length(Is);

    Ui = [];
    Ti = []; % Tracer
%
% Interpolate onto Z
    for j=1:nns    % individ segments
						dh = abs(dHn(:,j));
						dh(isnan(dh))=0;
						Ib=min(find(abs(dh)<1e-3));
						if Ib==1 | isnan(Un(1,j)),  % V or U on C grid near coast can be 0, p-point is not 0
								ti=ZZf*nan;
								si=ZZf*nan;
								ui=ZZf*nan;
						else
								zz=-cumsum(dh);
								t=Tr(:,j);
								u=Un(:,j);
								hb=-sum(dh);
								ibz = max(find(ZZf>=hb));
								nl=length(zz);
								for kl=Ib:nl
										zz(kl)=zz(kl-1)-0.1;
										t(kl)=t(kl-1);
										u(kl)=u(kl-1);
								end;
								zz=[0;zz];
								t=[t(1);t];
								u=[u(1);u];

								if abs(zz(end))<abs(ZZf(end))
										zz(end)=ZZf(end);
								end

								Ip=find(isnan(t));
								if ~isempty(Ip),
										fprintf('Interpolation fields: nans\n');
										keyboard
								end

								ti = interp1(zz,t,ZZf,'pchip');
								ti(ibz+1:end)=nan;
								ui = interp1(zz,u,ZZf,'pchip');
								ui(ibz+1:end)=nan;
						end
						Ti(:,j)=ti;
						Ui(:,j)=ui;
   	end % segments

    SCT(isct).Time(cc,1) = dnmb;
    SCT(isct).Unrm(cc,:,:)     = Ui;
    SCT(isct).Trcr(cc,:,:)     = Ti;
    

  end;  % sections

		fprintf('++++++>  Processed 1 record %8.5f min\n\n',toc/60);

		if s_mat>0 & mod(cc,5)==0
				fprintf('Saving %s\n',fmatout);
				save(fmatout,'SCT');
		end

end  % time

if s_mat>0 & cc>Nlast
		fprintf('Saving %s\n',fmatout);
		save(fmatout,'SCT');
end

