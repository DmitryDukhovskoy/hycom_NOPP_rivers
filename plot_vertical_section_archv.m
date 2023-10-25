% Vertical sections output archive files
% either nest files (instanteneous archv)
% or model output (saved as daily average archm)
addpath /usr/people/ddmitry/codes/MyMatlab/;
addpath /usr/people/ddmitry/codes/MyMatlab/hycom_utils;
addpath /usr/people/ddmitry/codes/MyMatlab/colormaps;
addpath /usr/people/ddmitry/codes/anls_mtlb_utils/hycom_arc08;
startup;

close all
clear

f_fort = 1; % =0 - check interp. fields created remap32to41lrs_fixedZ.m
            % =1 - created in remap_32to41FORT/remap32to41lrs.F90 

s_fig = 0;
%pfld  = 'temp';
%pfld  = 'salin';
pfld = 'v-vel.';
%pfld = 'u-vel.';
ptot = 0; % =1 - plot total; =0 - baroclinic u-vel/v-vel in archv

% GOFS3.0 vs GOFS3.1 nest fields:
gofs = 'gofs31'; % gofs30: for nest from gofs3.0 (32L, T07) and 
                 % gofs31: gofs3.1 (41L, T11)

R = 'ARCc0.08';
%R = 'ARCc0.04'; <--- Use plot_nest004_vertical_section.m

hg    = 2^100;
rg    = 9806;
% Change Tv and nlev for different ARCc fields
%Tv    = 7; % bathym version: 07 or 11
%nlev  = 32; % 32 or 41
%Tv    = 11; % bathym version: 07 or 11
%Tv    = 17; % bathym version: 07 or 11
nlev  = 41; % 32 or 41

expt = 110;  % ARCc HYCOM experiment #

yr=2004;
iday=362;
%iday=206; % 2005, 206 - fast-corrected on Gordon
hr=0;

j0=2134;
i0=896;

%sct = 'WE_Arct';
%sct = 'Atl_OB';
sct = 'pacif_ob';
%sct = 'pacif_merid';
switch(lower(sct))
 case('we_arct');
% West-East Central Arctic:
  jj1=1430;
  jj2=jj1;
  ii1=200;
  ii2=1450;
 case('atl_ob')
% Atl. OB:
  jj1=4;
  jj2=jj1;
%  ii1=190;
  ii1=1;
  ii2=1050;
 case('pacif_ob');
% Pacific OB:
  jj1=2519;
  jj2=jj1;
  ii1=1;
  ii2=1600;
 case ('pacif_merid');
  jj1=1524;
  jj2=2520;
  ii1=644;
  ii2=644;
end

if ii1==ii2
  xsct = 1; % section along X axis
else
  xsct = 0;
end


pthbin  = '/Net/kronos/ddmitry/hycom/ARCc0.08/tmp_files/';
pthrlx41= '/Net/kronos/ddmitry/hycom/ARCc0.08/relax_41layers/output/';
%ptharc41= '/Net/kronos/ddmitry/hycom/ARCc0.08/archv_41layers/';
ptharc41= sprintf('/Net/mars/ddmitry/hycom/%s/nest_files/',R);
%pthtopo = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/topo_grid/';
pthtopo = sprintf('/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/%s/topo_grid/',R);
pthfig  = '/nexsan/people/ddmitry/Net_ocean/HYCOM/ARCc/ARCc0.08/110/fig_vert_grid/';

switch(R)
 case('ARCc0.08'),
  Tv  = 11; % arc08
  TVn = '11';
  flnmbs = 'archv_arcT11L41';
 case('ARCc0.04'),
  Tv  = 17; % arc04
  TVn = '17DD';
%  flnmbs = 'archv_arc04T17L41'
  flnmbs = 'archv_arc04T17L41'
end

% Topo:
fltopo=sprintf('%sdepth_%s_%s.nc',pthtopo,R,TVn);
%fltopo=sprintf('%sdepth_ARCc0.04_%2.2iDD.nc',pthtopo,Tv);
HH  = nc_varget(fltopo,'Bathymetry');
LAT = nc_varget(fltopo,'Latitude');
LON = nc_varget(fltopo,'Longitude');
[mm,nn] = size(HH);
%HH(HH>0)=nan;


%iday=54;
if nlev == 41
% Check dummy relax files, there is no U,V in it
%  fina = sprintf('%sarchv.rlx.ARCc0.08%2.2i_%4.4i_%3.3i.a',pthrlx41,Tv,yr,iday);
%  finb = sprintf('%sarchv.rlx.ARCc0.08%2.2i_%4.4i_%3.3i.b',pthrlx41,Tv,yr,iday);
% Interpolated archive fields to 41 layers in matlab:
  fina = sprintf('%sarchv_arcT%2.2iL41.%i_%3.3i_00.a',ptharc41,Tv,yr,iday);
  finb = sprintf('%sarchv_arcT%2.2iL41.%i_%3.3i_00.b',ptharc41,Tv,yr,iday);

  if f_fort==1 % interpolated in Fortran:
% Diagnostic of U/V errors in Pacific sector in the nest files:
% *FORTerr - fields with errors, U/V are not rotated, i.e. 
%            opposite to what it should be in ARCc
% *FORT    - fields obatined after the bug is fixed in the 1st 
%            step when GLBb fields are rotated and remapped into ARCc
% *crct    - quick fix of the nest: U/V in the error nest fields are simply
%            witched by 180 degrees in the upper part of the ARCc grid
%            does not see any big difference from the "right" way
% ARCc0.04 files - prepared on newton
%                  from GOFS3.1 nest files ARCc0.08
%                  nan's (-999) at the northern OB are corrected
%                  (files with "fxd" at the end)
%
%   fina = sprintf('%sarchv_arcT%2.2iL41.%4.4i_%3.3i_00_crct.a',... % corrected fast
%		 ptharc41,Tv,yr,iday);
%    finb = sprintf('%sarchv_arcT%2.2iL41.%4.4i_%3.3i_00_crct.b',...
%		 ptharc41,Tv,yr,iday);
% fina = sprintf('%sarchv_arcT%2.2iL41.%4.4i_%3.3i_00.a',... % nest files from Newton
%		 ptharc41,Tv,yr,iday);
% finb = sprintf('%sarchv_arcT%2.2iL41.%4.4i_%3.3i_00.b',...
%		 ptharc41,Tv,yr,iday);
if strcmp(R,'ARCc0.04')
 fina = sprintf('%sarchv_arc04T%2.2iL41.%4.4i_%3.3i_00.a',... % nest files from Newton
		 ptharc41,Tv,yr,iday);
 finb = sprintf('%sarchv_arc04T%2.2iL41.%4.4i_%3.3i_00.b',...
		 ptharc41,Tv,yr,iday);
else
 fina = sprintf('%sarchv_arcT%2.2iL41.%4.4i_%3.3i_00.a',... % nest files from Newton
		 ptharc41,Tv,yr,iday);
 finb = sprintf('%sarchv_arcT%2.2iL41.%4.4i_%3.3i_00.b',...
		 ptharc41,Tv,yr,iday);
end  
% fina = sprintf('%sarchv_arcT%2.2iL41.%4.4i_%3.3i_%s.a',... % nest files from Newton
%		 ptharc41,Tv,yr,iday,gofs);
% finb = sprintf('%sarchv_arcT%2.2iL41.%4.4i_%3.3i_%s.b',...
%		 ptharc41,Tv,yr,iday,gofs);
%   fina = sprintf('%sarchv_arcT%2.2iL41.%4.4i_%3.3i_00_FORT.a',... % corrected right
%		 ptharc41,Tv,yr,iday);
%    finb = sprintf('%sarchv_arcT%2.2iL41.%4.4i_%3.3i_00_FORT.b',...
%		 ptharc41,Tv,yr,iday);
%    fina = sprintf('%sarchv_arcT%2.2iL41.%4.4i_%3.3i_00_FORTerr.a',...
%		 ptharc41,Tv,yr,iday);
%    finb = sprintf('%sarchv_arcT%2.2iL41.%4.4i_%3.3i_00_FORTerr.b',...
%		 ptharc41,Tv,yr,iday);
  end
  
elseif nlev == 32
% Check ARCc Topo11 or Topo 07 
% interpolated from GLBbT07, both 32 layers
  fina = sprintf('%sarchv_arcT%2.2i.%4.4i_%3.3i_%2.2i.a',pthbin,Tv,yr,iday,hr);
  finb = sprintf('%sarchv_arcT%2.2i.%4.4i_%3.3i_%2.2i.b',pthbin,Tv,yr,iday,hr);
end

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
%keyboard

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


for ik=1:2
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

  stt{1}=sprintf('%s, %s, ARCc0.08T%2.2iL%i from GLBb0.08T07L32, %i/%3.3i',...
	     sct,pfld,Tv,l,yr,iday);
  if ~isempty(Abtrop)
    stt{1}=sprintf('%s, %s+%s, ARCc0.08T%2.2iL%i from GLBb0.08T07L32, %i/%3.3i',...
	     sct,pfld,fadd,Tv,l,yr,iday);
  end  
  stt{2}=sprintf('%s',fina);
%  title(stt,'Interpreter','none');
  htt = text(0.1,10,stt,'Fontsize',13,'Interpreter','none');
  set(htt,'Position',[XL(1,1) 250]);

  txtbt='plot_vertical_section_archv.m';
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

 
