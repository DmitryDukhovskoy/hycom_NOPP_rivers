% Extract T, S, normal U at the specified straits
% Use daily mean output archm  - Total velocities !!!
% no need to add Ubarotropic
%
% Save 2D arrays (depth x width) every N days
%  Corrected flux calculation with T,S,U allocation is used
addpath /usr/people/ddmitry/codes/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/MyMatlab/Wavelet;
startup;

close all
clear

%TV=11;
YR1=2017;
YR2=2017;
dday=5;

ixx    = 10; % experiment name and dir - check with EXPT - expt 023
%ixx    = 6;  % expt 022 original 
EXPT   = sub_cice_experiments;
expt   = EXPT(ixx).Nmb;
texpt  = EXPT(ixx).cice_opt; % CICE options for sens. experiments
res    = EXPT(ixx).res;

if res == 0.04;
  fprintf(' Use extr_TSVdailymean_straits04.m for 0.04\n');
  error(' Res 0.04 not supported');
end
 
% CHeck if archm - do not add U/V barotropic !!!
s_mat=1;  % =2 - load saved and start from the last record

Cp    = 4200; % J/kg K
Tref1 = -1.8; % Ref T to calc. H flux
Tref2 = 0;    % Ref T to calc. H flux
Sref1 = 34.8;
Sref2 = 34.9; 
hgg   = 1e20;
f_zgrd=0;  % =1 - calculate fluxes from z-grid interpolated U,T,S - less accurate
           % mostly for comparison and validation

pthout = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_straits/';
pthmat = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/strait_fluxes/',expt);
pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';

btx='extr_TSVdailymean_straits08.m';

fprintf('arc08-%3.3i Heat and Vol fluxes Straits %i-%i, save=%i\n',...
        expt,YR1,YR2,s_mat);

ftopo = sprintf('%sdepth_ARCc0.08_11.nc',pthtopo); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);

% Fram Section is close to Moorings ~79N
%SCT = sub_define_sections04(HH,LON,LAT);
SCT = sub_define_AO_NA_sections08(HH,LON,LAT);
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
     1 0.4 0.6; ...
     0.45 0.2 0.85];


f_map=0;
if f_map==1
  fprintf('Drawing map with segments\n');
  fn=1;
%  sub_plot_Greenl_contour(HH,LON,LAT,fn,GC);
  figure(1); clf;
  contour(HH,[0 0],'k');
  hold on;
  contour(HH,[-5000:500:-100],'Color',[0.9 0.9 0.9]);

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
  set(gca,'xlim',[600 2800],...
          'ylim',[200 2800]);

  bottom_text(btx,'pwd',1);

  keyboard
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

  fmatout=sprintf('%shycom008_%3.3i_%s_StraitFluxesDayMean_%4.4i.mat',...
                    pthmat,expt,texpt,YR);

  cc=0;

% Daily Heat, and vol fluxes
  VHFLX=struct;
  nsct=length(SCT);
		for isct=1:nsct
				VHFLX(isct).Tref1=Tref1;
				VHFLX(isct).Tref2=Tref2;
				VHFLX(isct).Sref1=Sref1;
				VHFLX(isct).Sref2=Sref2;
				VHFLX(isct).Vfm=0;
				VHFLX(isct).Hf1m=0;
				VHFLX(isct).Hf2m=0;
				VHFLX(isct).FWf1=0;
				VHFLX(isct).FWf2=0;
				VHFLX(isct).T=[];
				VHFLX(isct).V=[];
				VHFLX(isct).U=[];
		end

  dlast=0;
  if s_mat == 2 
    fprintf('Loading saved %s start from last recrod\n',fmatout);
    load(fmatout);
    dlast = SCT(end).Time(end);
    fprintf('Last saved record %s\n\n',datestr(dlast));
%    keyboard
  end

% Find 1st record that exists in the directory:
  iday1 = 1;
  for iday=1:ndays
    dnmb=datenum(dJ1)+iday-1;
%    iday = dnmb-datenum(yr,1,1)+1;
    pthbin=sprintf('/nexsan/people/ddmitry/hycom/ARCc0.08_%3.3i/data/%4.4i_mean_%s/',...
                    expt,yr,texpt);
    fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
    finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);

    if exist(fina,'file') 
      fprintf('First file exist in dir, day=%i\n',iday);
      iday1=iday;
      break;
    end
  end    

  for iday0=iday1:dday:ndays
%  for iday=282:282
    dnmb=datenum(dJ1)+iday0-1;
    DV   = datevec(dnmb);
%    iday = dnmb-datenum(yr,1,1)+1;
    iday = iday0;

%    pthbin=sprintf('/nexsan/hycom/ARCc0.04_%3.3i/data/%4.4i/',expt,yr);
    pthbin=sprintf('/nexsan/people/ddmitry/hycom/ARCc0.08_%3.3i/data/%4.4i_mean_%s/',...
										expt,yr,texpt);

    fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
    finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);

    if ~exist(fina,'file') | ~exist(finb,'file')
      fprintf('Not found %s or %s\n\n',fina,finb);
      if dday==1; continue; end;
%
% Search for close days:
      for ik=1:dday
        iday = iday-1;
        dnmb = datenum(dJ1)+iday-1;
        DV   = datevec(dnmb);
        fina = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.a',pthbin,expt,yr,iday);
        finb = sprintf('%s%3.3i_archm.%4.4i_%3.3i_12.b',pthbin,expt,yr,iday);
        if exist(fina,'file'); 
          fprintf(' Found closest file: %s\n',fina);
          break; 
        end;
      end
% If went back to iday0 - skip
      if iday==iday0;
        continue;
      end
    end

    cc=cc+1;

    if dnmb<=dlast
      fprintf('Record exist %s, skipping ...\n',datestr(dnmb));
      continue
    end

    fprintf('Reading %4.4i/%2.2i/%2.2i: %s\n',DV(1:3),fina);

    tic;
    [F,n,m,nlr] = read_hycom(fina,finb,'temp');
    F(F>hgg)=nan;
    T=F;

    [F,n,m,nlr] = read_hycom(fina,finb,'salin');
    F(F>hgg)=nan;
    S=F;
%
    [F,n,m,nlr] = read_hycom(fina,finb,'u-vel.');
    F(F>hgg)=0;
    U=F;

    [F,n,m,llr] = read_hycom(fina,finb,'u_btrop');
    F(F>hgg)=0;
    Ub=squeeze(F);
%
%    for ilr=1:nlr
%      U(ilr,:,:)=squeeze(U(ilr,:,:))+Ub;
%    end

% These are total velocities    - archm
% barotropic NOT needed
    [F,n,m,nlr] = read_hycom(fina,finb,'v-vel.');
    F(F>hgg)=0;
    V=F;

    clear F

    [F,n,m,llr] = read_hycom(fina,finb,'v_btrop');
    F(F>hgg)=0;
    Vb=squeeze(F);
%
%    for ilr=1:nlr
%      V(ilr,:,:)=squeeze(V(ilr,:,:))+Vb;
%    end

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
       Si = [];
       Ti = [];
       UTS = struct;
       VFlx_btrp = [];
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
% Collocate - not collocated U,V,T
% clc = check_collocatedU(HH,U,V);
         if y1==y2  % horiz segment, V flux
           CLC=sub_collocate2uv(V,T,S,dH,HH,i0,j0,i1,j1);
         else
           CLC=sub_collocate2uv(U,T,S,dH,HH,i0,j0,i1,j1);
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
% For checking: calculate barotropic transport: 
         if y1==y2
           Unb = Vb(j0,i0);
         else
           Unb = Ub(j0,i0);
         end
         VFlx_btrp(j) = snrm*Unb*CLC.Hn*segm_dL(j);

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
% Plot layers:
          hold on
          for ilr=1:nlr
            zp0=zzp(ilr,:);
            plot(xxl,zzp,'k-');
          end
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
      hf1=nansum(un.*Cp.*rhow.*(tn-Tref1).*DX.*DZ);
      Hflx1=nansum(hf1);
      hf2=nansum(un.*Cp.*rhow.*(tn-Tref2).*DX.*DZ);
      Hflx2=nansum(hf2);
      fwf1=nansum(un.*(Sref1-sn)./Sref1.*DX.*DZ);
      FWflx1=nansum(fwf1);
      fwf2=nansum(un.*(Sref2-sn)./Sref2.*DX.*DZ);
      FWflx2=nansum(fwf2);
      I=find(sn<=Sref1);
      FWflx10=nansum(un(I).*(Sref1-sn(I))./Sref1.*DX(I).*DZ(I));

% Update:
      SCT(isct).Time(cc,1)       = dnmb;
      SCT(isct).VolFlx_m3s(cc,1) = Vflx;
      SCT(isct).Hflx1_W(cc,1)    = Hflx1;
      SCT(isct).Hflx2_W(cc,1)    = Hflx2;
      SCT(isct).FWflx1_m3s(cc,1) = FWflx1;
      SCT(isct).FWflx2_m3s(cc,1) = FWflx2;
%      SCT(isct).FWflx1_s1_m3s(cc,1) = FWflx10; % FW flux integrated only where S<Sref
      SCT(isct).Unrm(cc,:,:)     = Ui;
      SCT(isct).T(cc,:,:)        = Ti;
      SCT(isct).S(cc,:,:)        = Si;

% ==========================================
%  Diagnostics
% ==========================================
      Vfm  = mean(SCT(isct).VolFlx_m3s);
      hf1m = mean(SCT(isct).Hflx1_W);
      hf2m = mean(SCT(isct).Hflx2_W);
      fw1m = mean(SCT(isct).FWflx1_m3s);
      fw2m = mean(SCT(isct).FWflx2_m3s);
%      fw3m = mean(SCT(isct).FWflx1_s1_m3s);

      VFb = [];
      if exist('VFlx_btrp','var')
        VFb = nansum(VFlx_btrp);
      end
      
      fprintf('*****    \n');
      fprintf('Section %s\n',SCT(isct).Name);
      fprintf('Mean fluxes: Vol=%5.1f Sv, Heat1=%6.2f TW, Heat2=%6.2f TW\n',...
              Vfm*1e-6, hf1m*1e-12, hf2m*1e-12);
      if ~isempty(VFb)
        fprintf('Barotropic Vol Flux: %5.1f Sv\n',VFb*1e-6);
      end
      fprintf('FWFlux(34.8)=%5.1f mSv,  FWF(34.9)=%5.1f mSv\n',...
              fw1m*1e-3, fw2m*1e-3);

%keyboard
%      if isct==3; keyboard; end;

    end    % sections
    fprintf('===================\n');
    fprintf('++++++>  Processed 1 record %8.5f min\n\n',toc/60);


    if s_mat>=1 & mod(cc,2)==0
      fprintf('Saving %s\n',fmatout);
      save(fmatout,'SCT');
    end


  end  % day

end    % year

if s_mat>=1
		fprintf('Saving %s\n',fmatout);
		save(fmatout,'SCT');
end

