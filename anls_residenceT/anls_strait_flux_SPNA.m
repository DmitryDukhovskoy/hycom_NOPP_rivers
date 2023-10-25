% To test hypothesis about convergence/divergence
% of particles in the SPNA at different depths
%
% Estimate divergence - volume averaged
% using Gauss (diverg) theorem: integrate fluxes over the 
% boundaries
% Time series extracted in calc_flux_SPNA.m
%
addpath /Net/Movies0/ddmitry/MyMatlab
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /home/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers;
%addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
%addpath /Net/Movies0/ddmitry/MyMatlab/seawater
startup

close all
clear


pthmat = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/110/data_GG_prt/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.08/topo_grid/';
btx = 'anls_strait_flux_SPNA.m';

SLR = [10; 15; 23; 31]; % layer where Lagr. Part are advected
nlr = length(SLR);


YR1=2000;
YR2=2010;

ftopo = sprintf('%sdepth_ARCc0.08_11.nc',pthtopo); % 
fprintf('Getting topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);


% Get sections:
fsct = sprintf('%sSCT_SPNA.mat',pthmat);
fprintf('Loading %s\n',fsct);
load(fsct);


nsct = length(SCT);
f_map=0;
if f_map==1
  LMSK = HH*0;
  LMSK(HH<0)=1;
  cmph = [0 0 0; 1 1 1];

  figure(10); clf;
%  contour(HH,[0 0],'k');
  pcolor(LMSK); shading flat;
  caxis([0 1]);
  colormap(cmph);

  hold;
  contour(HH,[-1000 -1000],'Color',[0.5 0.5 0.5]);
  contour(HH,[-500 -500],'Color',[0.5 0.5 0.5]);

  for isg=1:nsct
    II = SCT(isg).I;
    JJ = SCT(isg).J;

    iocn = SCT(isg).Nocn;
    if iocn>0
      clr = [0.8 0.8 0.8];
    else
      clr = [0.8 0.8 0.8];
    end

    II = SCT(isg).I;
    JJ = SCT(isg).J;

    plot(II,JJ,'.-','Color',clr);
  end
  axis('equal');
  set(gca,'xlim',[300 1150],...
            'ylim',[50 800]);

  title('SPNA region');
  bottom_text(btx,'pwd',1);
end

% Pick ocean sections:
clear VFLX
iosc = 0;
jtot = 0;
for isc=1:nsct
  II = SCT(isc).I;
  JJ = SCT(isc).J;
  II = II(:);
  JJ = JJ(:);

		iocn = SCT(isc).Nocn;
		if iocn>0
    iosc=iosc+1;
    VFLX(iosc).II=II;
    VFLX(iosc).JJ=JJ;
    VFLX(iosc).Hb=SCT(isc).Hb;
    VFLX(iosc).Istart=jtot+1;
    VFLX(iosc).Iend=jtot+length(II);
  

    if f_map>0
      clr=rand(1,3);
      plot(II,JJ,'.-','Color',clr);
      text(mean(II),mean(JJ),sprintf('%i',iosc),'Fontsize',14);
    end
		end

  jtot=jtot+length(II)-1;
end 

%for io=1:length(VFLX);
%  i1=VFLX(io).Istart;
%  i2=VFLX(io).Iend;
%
%  plot([i1 i2],[0 0],'r-');
%end

% Mean fluxes by layers
% Calculate fluxes through open segements at different depth levels
tcc=0;
nocn = length(VFLX);
TM=[];
VFlx=[];
for YR=YR1:YR2
  fmat = sprintf('%shycom008_112_SPNAcntr_VolFlux_%i.mat',pthmat,YR);
  fprintf('Loading %s\n',fmat);
  load(fmat);

  tmm = GVHFLX.Time';
% Overall flux: checking
  fmm=GVHFLX.Vflx'; % vol transport total, daily m/3s depth-integrated
  VFlx=[VFlx;fmm];

  for icc=1:nocn
    i1=VFLX(icc).Istart;
    i2=VFLX(icc).Iend;
    A  = squeeze(GVHFLX.VflxPnts(:,:,i1:i2));
%    DZ = squeeze(GVHFLX.DZ(:,:,i1:i2)); 
% Total flux and by layers:
% + is into the box
    nt=length(tmm);
    tcc=length(TM);
    for it=1:nt
      tcc=tcc+1;
      dmm=squeeze(A(it,:,:));
      VFLX(icc).VFtot(tcc,1)=nansum(nansum(dmm)); % Vol Transp depth-intgr by sections
      Iout=find(dmm<0);  % outflow points
      VFLX(icc).Vout(tcc,1)=nansum(dmm(Iout));
      dmm_neg = dmm;
      dmm_neg(dmm>0)=0;
      for ilr=1:nlr
        if ilr>1
          il1=SLR(ilr-1)+1;
          il2=SLR(ilr);
        else
          il1=1;
          il2=SLR(ilr);
        end
        VFLX(icc).VFlr(ilr,tcc)=nansum(nansum(dmm(il1:il2,:)));  % transport in section, by layers
        VFLX(icc).VFlr_out(ilr,tcc)=nansum(nansum(dmm_neg(il1:il2,:))); 
      end
    end
  end
  TM = [TM;tmm];

end
% Approximate depth of last layer with Lagr floats, ilr=31
for icc=1:nocn
  i1=VFLX(icc).Istart;
  i2=VFLX(icc).Iend;
  DZ = squeeze(GVHFLX.DZ(:,:,i1:i2)); 
  tm0 = size(DZ,1);
  zL0 = [];
  ilr = SLR(end);
  for jtm=1:tm0
    amm = squeeze(DZ(jtm,:,:));
    zL  = -cumsum(amm,1);
    jz  = find(zL(ilr,:)==min(zL(ilr,:)),1);
    zL0(jtm) = 0.5*(zL(ilr,jz)+zL(ilr+1,jz)); % bottom interface
  end
  VFLX(icc).Z_LastLagrLr    = mean(zL0); % mean deepest value
  VFLX(icc).Z_DeepestBottom = min(VFLX(icc).Hb);
end

  aa = squeeze(DZ(:,ilr,:));


%
% Check total flux:
% sum by sections - daily depth-ingtgr flux should = total daily depth intgr 
nrc=length(TM);
Ftot=zeros(nrc,1);
for isc=1:nocn
 Ftot=Ftot+VFLX(isc).VFtot;  % daily mean vol flux m3/s
end;

fprintf('Long-term mean overall tranport %8.5f Sv, Ftot=%8.5f Sv\n',mean(VFlx*1e-6),mean(Ftot*1e-6));

% Report mean Vol fluxes by segments:
fprintf(' Mean vol transport over %i-%i:  \n',YR1,YR2);
for isc=1:nocn
  dmm = VFLX(isc).VFtot;
  nmm = VFLX(isc).Vout;
  fprintf('----------------------------------\n');
  fprintf('Sect %i, Depth-intgr Overall: Net  %12.8f Sv/ outfl: %12.8f Sv \n',...
           isc,mean(dmm)*1e-6, mean(nmm)*1e-6);
  fprintf(' Min Z(lr= %i) = %8.2f m , Deepest pnt %8.2f m\n', ...
           SLR(end),VFLX(isc).Z_LastLagrLr, VFLX(isc).Z_DeepestBottom);

  lrtot=0;
  lrout=0;
  for ilr=1:nlr
				if ilr>1
						il1=SLR(ilr-1)+1;
						il2=SLR(ilr);
				else
						il1=1;
						il2=SLR(ilr);
				end

    dmm = VFLX(isc).VFlr(ilr,:);
    nmm = VFLX(isc).VFlr_out(ilr,:);
    fprintf('      Layers %i-%i   : Net %12.8f Sv/ outfl: %12.8f Sv \n',...
              il1,il2,mean(dmm)*1e-6, mean(nmm)*1e-6);
    lrtot=lrtot+mean(dmm)*1e-6;
    lrout=lrout+mean(nmm)*1e-6;
  end
  fprintf('      All Layers 1-31 : Net %12.8f Sv/ outfl: %12.8f Sv \n',lrtot,lrout); 
end

  



