% Plot vertical disctribution
% of Greenland tracer
% on South Gr shelves
% Sections extracted in xsection_tracer_sGrSh008.m
%
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/seawater;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

nTr   = 1;
YR1 = 2016;
YR2 = 2016; 

plr=0; % highlight this interface
btx = 'plot_xsection_tracer_fjords008.m';


regn = 'ARCc0.08';
%expt = 110; % experiment without runoff
expt = 112;  % epxeriment with Greenland runoff and monthly Arctic rivers
pthmat  = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt);
pthtopo= '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthriv = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/110/data_mat/';

ftopo = sprintf('%s/depth_%s_11.nc',pthtopo,regn); % 
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');
[mm,nn]=size(LON);
[DX,DY]=sub_dx_dy(LON,LAT);


SGM  = sub_sGrShelf_sections(HH,LON,LAT);
nsct = length(SGM);


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
  set(gca,'xlim',[400 960],...
	  'ylim',[350 1070]);

  btx = 'plot_xsection_tracer_sGrSh008.m';
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


fmatout = sprintf('%shycom008_%3.3i_Trcr%2.2i_sGrShxsct_%4.4i-%4.4i.mat',...
                    pthmat,expt,nTr,YR1,YR2);
fprintf('Loading %s\n',fmatout);
load(fmatout);
nsct=length(SCT);


%
% COnvert Tr concentration to GFWA m3/m3 of sea water

frv = sprintf('%sGreenland_cumFWFlux.mat',pthriv);

f_griv = 0; %=1 - rederive cumulative Greenland FWFlux from Bamber
if f_griv==1
		sub_get_GrRunoff(friv)
end
fprintf('f_griv %i, Loading %s\n',f_griv,frv);
load(frv);

TM = SCT(1).Time;
imo0=-1;
for it=1:length(TM)
		DV = datevec(TM(it));
		iyr = DV(1);
		imo = DV(2);
  YR = iyr;
		if imo~=imo0
				expt0=110;
				if YR<1999
						pthmat2=sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/%3.3i/data_mat/',expt0);
				else
				pthmat2=sprintf('/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/%3.3i/data_matTr/',expt);
				end
				fmat2 = sprintf('%sMassTr%2.2i_lrs_%i%2.2i.mat',pthmat2,nTr,iyr,imo);
				fprintf('Loading %s\n',fmat2);
				load(fmat2);
				imo0=imo;
		end
% find whole-depth layer
		nlr = length(TRCR);
		ibtm=5; % whole-depth tracer mass
		Tr_dom=squeeze(TRCR(ibtm).MassTr_kg); % tracer integrated ove whole water depth
		amm=Tr_dom./(abs(HH).*DX.*DY);  % kg/m3
		Tr_dom(amm<=0.125)=nan;
		MTr_dom(it,1) = nansum(nansum(Tr_dom)); % overall mass of the tracer in the domain
end


% Plot sections
% average over year
for isc=1:nsct
  Hb = SCT(isc).Hbottom;
  Un = squeeze(nanmean(SCT(isc).Unrm,1));
  Tr = squeeze(nanmean(SCT(isc).Trcr,1)); % kg/m3 of Tracer
  Xdst = SCT(isc).dist_m*1e-3;
  ZZ = SCT(isc).ZZ;
  xU = Xdst;
  Hb = SCT(isc).Hbottom;
  TM = SCT(isc).Time;
  Is = SCT(isc).IJ_indx(:,1);
  Js = SCT(isc).IJ_indx(:,2);

% Convert Tr 
% Annual mean values assumed
% 1 year
		dv = datevec(TM);  
		ism = find(Ygr==dv(1,1));
		if isempty(ism), ism=length(cFWF); end; % after 2016 - same 
		fwf0 = cFWF(ism); % km3

		MTr_dom_mn=mean(MTr_dom);
		rr=Tr./MTr_dom_mn;
		GrVol=rr*fwf0*1e9;      % m3 of GFWA in 1 m3 of ocean, 1e-3 m3 = 1 L
  mVTr = GrVol;

% fill nans for plotting:
  for ipp=1:length(Is)
    dmm=mVTr(:,ipp);
    i1=max(find(~isnan(dmm)));
    if ~isempty(i1)
      dmm(i1+1:end)=dmm(i1);
    end
    mVTr(:,ipp)=dmm;
  end
  Un(isnan(Un))=0;

  lmVTr = log(mVTr);
  mUU   = Un;


  xBtm=[Xdst(1); Xdst; Xdst(end)];
  yBtm=[-6000; Hb; -6000];
  c1 = -7;
  c2 = -4.5;
  yl1= 1.02*min(Hb);
  xl1= 0;
  xl2= max(Xdst);

  switch(isc)
   case(1)
    CTR.Upos=[0.02:0.02:0.6];
    CTR.Uneg=[-0.6:0.02:-0.001];
   case(2)
    CTR.Upos=[0.05:0.05:0.8];
    CTR.Uneg=[-0.5:0.05:-0.001];
   case(3)
    CTR.Upos=[0.05:0.05:0.8];
    CTR.Uneg=[-0.8:0.05:-0.001];
   case(4)
    CTR.Upos=[0.05:0.05:0.8];
    CTR.Uneg=[-0.8:0.05:-0.001];
  end

  sub_plot_TrSct(isc,Xdst,ZZ,lmVTr,xBtm,yBtm,xU,mUU,c1,c2,yl1,xl1,xl2,CTR);
 
  nm=SCT(isc).Name;
  if strncmp(nm,'WGr',3)
    set(gca,'xdir','reverse');
  end

  sttl=sprintf('%s, ln(GFWA) m3 in 1 m3, log, %i-%i',nm,YR1,YR2);
  title(sttl);

  btx = 'plot_xsection_tracer_sGrSh008.m';
  bottom_text(btx,'pwd',1);

end
  
  

