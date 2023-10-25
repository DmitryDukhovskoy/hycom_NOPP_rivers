% Extract tracer  at the specified straits
% Save 2D arrays (depth x width) every N days
%  Corrected flux calculation with T,S,U allocation is used
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/Wavelet;
startup;

close all
clear


expt=112;
TV=11;
YR1=2001;
YR2=2002;
dday=7;
nTr=1;  % extract 1 tracer at a time, 1 - Greenland

s_mat=1;  % ==2 - start from last saved rec # in YEAR

hgg=1e20;
f_zgrd=0;  % =1 - calculate fluxes from z-grid interpolated U,T,S - less accurate
           % mostly for comparison and validation
ptharch= '/nexsan/people/ddmitry/';
pthout = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_straits/';
pthmat = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_theresa/';
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';

btx='extr_Trcr_daily_straits08.m';

fprintf('arc08-%3.3i Tracer Fluxes Gates %i-%i, save=%i\n',...
        expt,YR1,YR2,s_mat);


ftopo = sprintf('%sdepth_ARCc0.08_%2.2i.nc',pthtopo,TV); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);

% Fram Section is close to Moorings ~79N
SCT = sub_define_SPNA_sections(HH,LON,LAT);
nsct = length(SCT);


CLR=[0 0.4 0.8; ...
     0.8 0.4 0; ...
     1 0.2 0; ...
     0 1 0; ...
     0.8 0 0.6; ...
     0.8 1 0; ...
     0.3 0.7 1; ...
     1 0.7 0.5; ...
     0.9 0.3 0.6; ...
     0.2 0.7 0.4; ...
     0.1 0.4 0.3; ...
     0.8 0.3 0.9; ...
     1 0.4 0.6];

f_map=0;
if f_map==1
  fprintf('Drawing map with segments\n');
  fn=1;
%  sub_plot_Greenl_contour(HH,LON,LAT,fn,GC);
  figure(1); clf;
  contour(HH,[0 0],'k');
  hold on;
  contour(HH,[-5000:500:-100],'Color',[0.6 0.6 0.6]);

  for ip=1:nsct
%    IJ=SCT(ip).IJ;
%    plot([IJ(1,1) IJ(2,1)],[IJ(1,2) IJ(2,2)],...
%        'Linewidth',2.5,'Color',[1. 0.6 0]);
   clr=CLR(ip,:);
    IIs=SCT(ip).I;
    JJs=SCT(ip).J;
    plot(IIs,JJs,'-',...
         'Linewidth',2.5,'Color',clr);
    Ip=SCT(ip).IJPR(1);
    Jp=SCT(ip).IJPR(2);

    plot(Ip,Jp,'.','Markersize',14,'Color',clr);
  end

  axis('equal');
  set(gca,'xlim',[300 1300],...
          'ylim',[100 1300]);

  bottom_text(btx,'pwd',1);

%keyboard
end


ZZf = [(0:-1:-10)';(-12:-2:-200)';(-205:-5:-500)';...
       (-510:-10:-2000)';(-2025:-25:-2500)';(-2600:-100:-5000)'];
kzz = length(ZZf);

dZf=abs(diff(ZZf));
ZMf = [];
for ik=1:kzz-1
  ZMf(ik,1)=ZZf(ik)+0.5*(ZZf(ik+1)-ZZf(ik));
end




% Find contour through UV pnts for the sections
for isct=1:nsct
  IIs=SCT(isct).I;
  JJs=SCT(isct).J;
		nin=SCT(isct).IPR;
		UVGR = sub_UVpnts_contour(IIs,JJs,nin,HH);
		dL=sub_segm_dL(UVGR.sgmX,UVGR.sgmY,LON,LAT);

		SCT(isct).Norm=UVGR.Norm;
		SCT(isct).sgmX=UVGR.sgmX;
		SCT(isct).sgmY=UVGR.sgmY;
		SCT(isct).IJ_indx=UVGR.gridIndx_IJ; % grid indices corresponding to UV segments
		SCT(isct).adjIndx_I1J1=UVGR.adjIndx_I1J1; % adjacent grid point index - for collocating dH to UV pnts
		SCT(isct).segm_dL=dL;

  SCT(isct).ZZintrp = ZZf;
end

%keyboard
f_box=0;
if f_box==1
% Plot all boxes and sections and norms
% for each box
  fprintf('Drawing map with boxes\n');
  nf=15;
  sub_check_sct(nf,HH,SCT);
end


for YR=YR1:YR2
  yr=YR;
  dE=datenum(yr,12,31);
  dJ1=datenum(yr,1,1);
  ndays=dE-dJ1+1;

  fmatout=sprintf('%shycom008_%3.3i_Trcr%2.2i_StraitFluxesDay_%4.4i.mat',...
                    pthmat,expt,nTr,YR);

  cc=0;

% Daily Tr fluxes
  VHFLX=struct;
  nsct=length(SCT);
		for isct=1:nsct
    VHFLX(isct).VolFlx_m3s=0;  % vol flux
				VHFLX(isct).TrFlx=0;
				VHFLX(isct).TrcVrt=[]; % 2D Section of the tracer field
		end

  Nlast=0;
  if s_mat==2
    clear SCT
    fprintf('\n\n !!!!!!  Continue from the last saved record in %s !!!\n\n',fmatout);
    load(fmatout);

    Nlast = size(SCT(1).TrFlx,1);
    fprintf('Last saved record = %i\n\n',Nlast);

  end

  for iday=1:dday:ndays
%  for iday=282:282
    dnmb=datenum(dJ1)+iday-1;
    DV   = datevec(dnmb);
    iday = dnmb-datenum(yr,1,1)+1;
    pthbin = sprintf('/nexsan/archive/ARCc0.08_%3.3i/data/%4.4i/',expt,yr);
    if expt==112,
      pthbin=sprintf('/nexsan/hycom/ARCc0.08_%3.3i/data/%4.4i/',expt,yr);
    end

    fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
    finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);

    if ~exist(fina,'file') | ~exist(finb,'file')
      fprintf('Not found %s or %s\n\n',fina,finb);
      continue;
    end

    cc=cc+1;
    if cc<=Nlast
      fprintf('Record exist, skipping %4.4i/%2.2i/%2.2i \n',DV(1:3));
      for isc=1:nsct
        SCT(isc).Time(cc,1)=dnmb;  % add missing time array
      end
      continue;
    end

    fprintf('Reading %4.4i/%2.2i/%2.2i: %s\n',DV(1:3),fina);

    tic;
    nTr=1; % Greenland tracer
    [F,n,m,nlr] = read_hycom(fina,finb,'tracer','r_tracer',nTr);
    F(F>hgg)=nan;
    F(F<0)=0;
    Ctr=F;  % kg/m3
%
    [F,n,m,nlr] = read_hycom(fina,finb,'u-vel.');
    F(F>hgg)=0;
    U=F;

% These are total velocities    
% barotropic not needed
    [F,n,m,nlr] = read_hycom(fina,finb,'v-vel.');
    F(F>hgg)=0;
    V=F;

    [ZMh,ZZh] = sub_zz_zm(fina,finb,HH);
    dH=abs(diff(ZZh,1));

%
% Calculate fluxes through sections: 
    for isct=1:nsct
       iGTs=SCT(isct).sgmX;   % i indices of the uv segments start/end points
       jGTs=SCT(isct).sgmY;   % j indices of the uv segments
       IJp=SCT(isct).IJ_indx;               % grid index corrsponding to the segment
       IJadj=SCT(isct).adjIndx_I1J1;        % corresponding adjacent grid cell for collocation
       Nrm=SCT(isct).Norm;                  % unit norm in the box
       segm_dL=SCT(isct).segm_dL;           % segment length, m
       nns=length(iGTs(:,1));
       Ui = [];
       Ti = []; % Tracer
       Si = []; % Not needed
       UTS = struct;
       for j=1:nns
         x1=iGTs(j,1);
         x2=iGTs(j,2);
         y1=jGTs(j,1);
         y2=jGTs(j,2);
         xnrm=Nrm(j,1);
         ynrm=Nrm(j,2);
%
% Grid cell and adjacent cell
         i0=IJp(j,1);
         j0=IJp(j,2);
         i1=IJadj(j,1);
         j1=IJadj(j,2);
% Collocate
% Note subroutine is designde for 2 tracer fields, T/S, thus provide 2 Tracer arrays
         if y1==y2  % horiz segment, V flux
           CLC=sub_collocate2uv(V,Ctr,Ctr,dH,HH,i0,j0,i1,j1);
         else
           CLC=sub_collocate2uv(U,Ctr,Ctr,dH,HH,i0,j0,i1,j1);
         end
         In=find(CLC.dHn<1e-3);
         CLC.Tn(In)=nan;
         CLC.Sn(In)=nan;
         CLC.Un(In)=nan;

         snrm=sign(xnrm+ynrm); % one is 0,  need only direction
         UTS.Hb(j)=-CLC.Hn;
         UTS.Unrm(:,j)=snrm*CLC.Un;
         UTS.Tnrm(:,j)=CLC.Tn;
         UTS.Snrm(:,j)=CLC.Sn;
         UTS.dH(:,j)=CLC.dHn;

% 
% Interpolate onto Z levels
         dh=-CLC.dHn;
         dh(isnan(dh))=0;
         Ib=min(find(abs(dh)<1e-3));
         if Ib==1, 
           ti=ZZf*nan;
           si=ZZf*nan;
           ui=ZZf*nan;
         else
           zz=cumsum(dh);
           t=CLC.Tn;
           s=CLC.Sn;
           u=CLC.Un;
           hb=sum(dh);
           ibz = max(find(ZZf>=hb));
           nl=length(zz);
           for kl=Ib:nl
             zz(kl)=zz(kl-1)-0.1;
             t(kl)=t(kl-1);
             s(kl)=s(kl-1);
             u(kl)=u(kl-1);
           end;
           zz=[0;zz];
           t=[t(1);t];
           s=[s(1);s];
           u=[u(1);u];

           if abs(zz(end))<abs(ZZf(end))
             zz(end)=ZZf(end);
           end

           Ip=find(isnan(s));
           if ~isempty(Ip),
             fprintf('Interpolation fields: nans\n');
             keyboard
           end

           si = interp1(zz,s,ZZf,'pchip');
           si(ibz+1:end)=nan;
           ti = interp1(zz,t,ZZf,'pchip');
           ti(ibz+1:end)=nan;
           ui = interp1(zz,u,ZZf,'pchip');
           ui(ibz+1:end)=nan;
         end
         Si(:,j)=si;
         Ti(:,j)=ti;
         Ui(:,j)=ui;
      end % segments

% =========================
         fchck=0;
         if fchck==1
           DZ=-abs(UTS.dH);
           zzp=cumsum(DZ);
%           un=UTS.Unrm;
%           tn=UTS.Tnrm;
           sn=UTS.Snrm;
           
           figure(2); clf;
           axes('Position',[0.08 0.6 0.8 0.32]);
           xxl=[1:nns];
           pcolor(xxl,ZZf,Si); shading flat;
           set(gca,'ylim',[min(min(zzp)) 0]);
           title('Interpolated');
           colorbar

          axes('Position',[0.08 0.08 0.8 0.32]);
          pcolor(xxl,zzp,sn); shading flat;
           set(gca,'ylim',[min(min(zzp)) 0]);

          
         end
% =========================

      un=UTS.Unrm;
      tn=UTS.Tnrm;
      sn=UTS.Snrm;
      DZ=abs(UTS.dH);
      dx=segm_dL;

      [DX,dmm]=meshgrid(dx,[1:nlr]);

      rhow=sw_dens0(sn,tn);

% Fluxes over whole sections    
      vf0=nansum(un.*DX.*DZ);
      Vflx=nansum(vf0);
      TrFlx=nansum(un.*tn.*DX.*DZ);

% Update:
      SCT(isct).Time(cc,1)       = dnmb;
      SCT(isct).VolFlx_m3s(cc,1) = Vflx;
      SCT(isct).TrFlx(cc,1)      = nansum(TrFlx);
      SCT(isct).Unrm(cc,:,:)     = Ui;
      SCT(isct).Trcr(cc,:,:)     = Ti;

% ==========================================
%  Diagnostics
% ==========================================
      Vfm  = mean(SCT(isct).VolFlx_m3s);
      Trm  = mean(SCT(isct).TrFlx);
      
      fprintf('*****    \n');
      fprintf('Section %s\n',SCT(isct).Name);
      fprintf('Mean Vol flux = %6.3f Sv, Tracer Flx = %6.2g \n',...
              Vfm*1e-6, Trm);

    end    % sections
    fprintf('===================\n');
    fprintf('++++++>  Processed 1 record %8.5f min\n\n',toc/60);


    if s_mat>1 & mod(cc,10)==0
      fprintf('Saving %s\n',fmatout);
      save(fmatout,'SCT');
    end


  end  % day

		if s_mat>1 & cc>Nlast
				fprintf('Saving %s\n',fmatout);
				save(fmatout,'SCT');
		end
end    % year


