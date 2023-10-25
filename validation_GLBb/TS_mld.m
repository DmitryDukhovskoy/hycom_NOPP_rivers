addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
startup;

close all
clear

s_mat  = 1; % =1 - extract and save mat file; 
            % =0 - do not save; 
	    % = 2 - load existing year and fill missing months

yr1 = 1993;
yr2 = 2000;


rg=9806;  % convert pressure to depth, m
dT = 0.3; % T change to calculate d(rho) for MLD - Kara et al., 2003 & 2000

fprintf(' =====   Calculating MLD DDukhovskoy Method  =====\n');
fprintf(' =====   Time: %i - %i             =====\n\n',yr1,yr2);

pthbin  = '/Net/mars/ddmitry/hycom/GLBb2ARCc0.08/mnth_mean/'; 
pthmat  = '/Net/mars/ddmitry/hycom/GLBb0.08/data_mld/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/'; 
pthfig  = '/Net/mars/ddmitry/hycom/GLBb0.08/fig_MLD/';


TV = '07';
ftopo = sprintf('%sdepth_ARCc0.08_%s.nc',pthtopo,TV); % 
HH  = nc_varget(ftopo,'Bathymetry'); 
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(HH); 

% Mask of the region of interest
% exclude Pacific Ocean and North.Atl.
Lmsk             = HH*0;
Lmsk(HH<-500)    = 1;
Lmsk(1935:end,:)   = 0;  
Lmsk(1:400,1120:end)= 0;
Lmsk(1:500,1200:end)= 0;
%Lmsk(1:20,:) = 0;  % 

Iocn=find(Lmsk==1);
Ilnd=find(Lmsk==0); 


tmm=[-1.8:0.1:8];
smm=[30:0.1:35.1];
[TP,SL]=meshgrid(tmm,smm);
RHO=sw_dens0(SL,TP)-1000;

cmp = colormap(jet(32));

% 
for iyr=yr1:yr2
  year=iyr;
  MILD = struct; 
  im1 = 1;
  im2 = 12;
  MNTH=[1:12];
  
  for im=MNTH
    dnmb = datenum(iyr,im,1);
    dnmb2 = dnmb+35;
    dv2 = datevec(dnmb2);
    dnmb2 = datenum(dv2(1),dv2(2),1);
    
    fprintf('%i/%2.2i\n',year,im);  
    if year<1995  
      EE=190;
    elseif year>=1995 & year<2013
      EE=191;  
      if year == 1995 & im<=7
	EE=190;
      end;
    else
      EE1 = sub_exptGLBb(dnmb);
      EE2= sub_exptGLBb(dnmb2); % need the whole month 
      EE = EE2;
    end
    E=sprintf('GLBb0.08_%3.3i',EE);
    
    fina = sprintf('%s%3.3i_archMN_GLBb2ARCc.%i_%2.2i.a',...
		   pthbin,EE,year,im);      
    finb = sprintf('%s%3.3i_archMN_GLBb2ARCc.%i_%2.2i.b',...
		   pthbin,EE,year,im);      
    
    if ~exist(fina,'file');
      fprintf('does not exist %s\n',fina);
      continue;
    end
    
%    cc=cc+1;
    fprintf('Reading %s\n',fina);
    
% Get layer thickness to be able
% to construct depth arrays of model layers
    fld='thknss';
    [F,n,m,l] = read_hycom(fina,finb,fld);  
    F(F>1e20)=nan;                          
    DP=F./rg;                               
    DP(DP<0.1)=nan; % 0-m layers, vanished  
%
% Interface depths and depths of the middle of the layers (m)
% NOTE: sign convention: depths are negative
% dP - Pa or m
    [ZZ,ZM] = sub_thck2dpth(DP); 

    fld='salin';
    [F,n,m,l] = read_hycom(fina,finb,fld);
    F(F>1e20)=nan;                        
    S=F;  
    S(isnan(DP))=nan;

    fld='temp';
    [F,n,m,l] = read_hycom(fina,finb,fld); 
    F(F>1e20)=nan;                         
    Temp=F; 
    Temp(isnan(DP))=nan;
    
    ii=800;
    jj=350;
    t=squeeze(Temp(:,jj,ii));
    s=squeeze(S(:,jj,ii));
    rho=sw_dens0(s,t);
    zm = squeeze(ZM(:,jj,ii));
    zz = squeeze(ZZ(:,jj,ii));
    [t,s] = sub_stabilize_rho(t,s,zm,zz);
    rho=sw_dens0(s,t);
    
    keyboard
    
    figure(3); clf;
    hold
    for kk=1:l-1
      clr=cmp(kk,:);
      plot(rho(kk),zm(kk),'.','Color',clr,'markersize',12);
    end
    
    
    figure(4); clf;
    hold on;
    for kk=1:l-1
      clr=cmp(kk,:);
      plot(s(kk),t(kk),'.','Color',clr,'markersize',12);
    end;
    contour(SL,TP,RHO,[27.4:0.05:28],'Color',[0.5 0.5 0.5]);
    contour(SL,TP,RHO,[27.0:0.5:28.5],'c');
    set(gca,'xlim',[34.8 35.1],...
	    'ylim',[2 7]);
    
    
  end
end