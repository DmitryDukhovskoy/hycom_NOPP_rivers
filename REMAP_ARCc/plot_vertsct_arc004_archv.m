% Vertical sections output archive files ARCc0.04
% either nest files (instanteneous archv)
% or model output (saved as daily average archm)
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

s_fig = 0;
%pfld  = 'temp';
%pfld  = 'salin';
pfld = 'v-vel.';
ptot = 0; % =1 - plot total; =0 - baroclinic u-vel/v-vel in archv

% GOFS3.0 vs GOFS3.1 nest fields:
gofs = 'gofs31'; % gofs30: for nest from gofs3.0 (32L, T07) and 
                 % gofs31: gofs3.1 (41L, T11)

hg    = 2^100;
rg    = 9806;
% Change Tv and nlev for different ARCc fields
%Tv    = 7; % bathym version: 07 or 11
%nlev  = 32; % 32 or 41
Tv    = 17; % bathym version: 07 or 11
nlev  = 41; % 32 or 41

expt = 110;  % ARCc HYCOM experiment #

yr=2007;
iday=001;
%iday=206; % 2005, 206 - fast-corrected on Gordon
hr=0;

%j0=2134;
%i0=896;

%sct = 'WE_Arct';
%sct = 'Atl_OB';
sct = 'pacif_ob';
%sct = 'pacif_merid';
switch(lower(sct))
 case('we_arct');
% West-East Central Arctic:
  jj1=2860;
  jj2=jj1;
  ii1=400;
  ii2=2900;
 case('atl_ob')
% Atl. OB:
  jj1=4;
  jj2=jj1;
%  ii1=190;
  ii1=1;
  ii2=2100;
 case('pacif_ob');
% Pacific OB:
  jj1=5038;
  jj2=jj1;
  ii1=1;
  ii2=3200;
 case ('pacif_merid');
  jj1=3048;
  jj2=5040;
  ii1=1288;
  ii2=1288;
end

if ii1==ii2
  xsct = 1; % section along X axis
else
  xsct = 0;
end


ptharc41= '/Net/mars/ddmitry/hycom/ARCc0.04/nest_files/';
pthtopo = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/topo_grid/';
pthfig  = '/Net/ocean/ddmitry/HYCOM/ARCc/ARCc0.04/010/fig_vert_grid/';

ftopo = sprintf('%sdepth_ARCc0.04_17DD.nc',pthtopo); % 
fprintf('Reading topo %s\n',ftopo);
HH  = nc_varget(ftopo,'Bathymetry');
LON = nc_varget(ftopo,'Longitude');
LAT = nc_varget(ftopo,'Latitude');

% Diagnostic of U/V errors in Pacific sector in the nest files:
%
% ARCc0.04 files - prepared on newton
%                  from GOFS3.1 nest files ARCc0.08
%                  nan's (-999) at the northern OB are corrected
%                  (files with "fxd" at the end)
%
fina = sprintf('%sarchv_arc04T%2.2iL41.%4.4i_%3.3ifxd.a',... % nest files from Newton
		 ptharc41,Tv,yr,iday);
finb = sprintf('%sarchv_arc04T%2.2iL41.%4.4i_%3.3ifxd.b',...
		 ptharc41,Tv,yr,iday);

fprintf('Opening %s\n',fina);
Abtrop = [];
if strncmp(pfld(2:5),'-vel',4) & ptot>0
  fadd = [pfld(1),'_btrop'];
  fprintf('Adding U/V btrop field %s to baroclinic\n',fadd);
  [F,n,m,l] = read_hycom(fina,finb,fadd);
  F = squeeze(F);
  F(F>1e20)=nan;
  if ~xsct
    Abtrop = squeeze(F(jj1,ii1:ii2));
  else
    Abtrop = squeeze(F(jj1:jj2,ii1));
  end
  
end


[F,n,m,l] = read_hycom(fina,finb,pfld);
F(F>1e20)=nan;
%B=squeeze(F(1,:,:));
if ~xsct
  A=squeeze(F(:,jj1,ii1:ii2)); 
else
  A=squeeze(F(:,jj1:jj2,ii1)); 
end

[a1,a2]=size(A);
if ~isempty(Abtrop)
  for k=1:a1
    A(k,:)=A(k,:)+Abtrop;
  end
end

%keyboard
% Add 0 surface for plotting:
A=[A(1,:);A];


fld='thknss';
[F,n,m,l] = read_hycom(fina,finb,fld);
%keyboard
F(F>1e10)=nan;
F(F<0.1)=nan;
F=F./rg;
if ~xsct
  Dsec=squeeze(F(:,jj1,ii1:ii2)); 
else
  Dsec=squeeze(F(:,jj1:jj2,ii1)); 
end  

Dsec(Dsec==0)=nan;

% Create Depth array of interface depths:
% Note these are BOTTOM interfaces 
% So Layer 1 is between interfaces 0m and ZZ(1)
clear ZZb
Dsec(isnan(Dsec))=0;
ZZb=zeros(1,a2);  % surface 0m
ZZb(2,:)=-Dsec(1,:);
for kk=2:l
  kzb=kk+1;
  ZZb(kzb,:)=ZZb(kzb-1,:)-Dsec(kk,:);
end

Nintf=size(ZZb,1);
for ii=1:a2
%  I=ii1+ii-1;
  I=ii;
  XL(1:Nintf,ii)=ones(Nintf,1)+I;
end


for ik=1:1
  figure(ik); clf;
  pcolor(XL,ZZb,A); shading flat
  colormap(jet);
  hold on;
  caxis([-0.2 0.2]);
  h=colorbar;
  set(h,'YColor',[0 0 0]);

  if ik==2
    for kk=1:l
      zz=ZZb(kk+1,:); % Plot layer bottom interface
      xx=XL(kk,:);
      if mod(kk,5)==0;
	plot(xx,zz,'k-','Color',[0.4 0.4 0.4]);
      else
	plot(xx,zz,'k--','Color',[0.7 0.7 0.7]);
      end    

    end
  end;

  set(gca,'color',[0 0 0]);

  stt{1}=sprintf('%s, %s, ARCc0.04T%2.2iL%i from nest ARCc0.08 (GOFS3.1), %i/%3.3i',...
	     sct,pfld,Tv,l,yr,iday);
  if ~isempty(Abtrop)
stt{1}=sprintf('%s, %s+%s, ARCc0.04T%2.2iL%i from nest ARCc0.08 (GOFS3.1), %i/%3.3i',...
	     sct,pfld,fadd,Tv,l,yr,iday);
  end  
  stt{2}=sprintf('%s',fina);
%  title(stt,'Interpreter','none');
  htt = text(0.1,10,stt,'Fontsize',13,'Interpreter','none');
  set(htt,'Position',[XL(1,1) 250]);

  txtbt='plot_vertsct_arc004_archv.m';
  bottom_text(txtbt,'Position',[0.06 0.06 0.8 0.05],'pwd',1);

  if s_fig>0
    fgnm=sprintf('%sARCcT%2.2i_L%i_%s_%s_%i_%3.3i-%2.2i',...
		 pthfig,Tv,nlev,pfld,sct,yr,iday,ik);
    if f_fort==1
      fgnm=sprintf('%sARCcT%2.2i_L%i_%s_%s_%i_%3.3i-%2.2iFORT',...
		 pthfig,Tv,nlev,pfld,sct,yr,iday,ik);
    end
    fprintf('Saving %s\n',fgnm);
    print('-dpng','-r150',fgnm);
  end

end;

 
