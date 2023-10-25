function SCT = sub_define_sections(HH,LON,LAT);
%
%
% Use Greenland "Boxes" - from heat fluxes Theresa
%pthmat = '/nexsan/people/ddmitry/Net_tholia/hycom/ARCc0.08/data_theresa/';
%YR=2005;
%fmatout=sprintf('%shycom008_%3.3i_Greenl_flx_POPboxes_%4.4i.mat',...
%                 pthmat,expt,YR);
%load(fmatout);



SCT = struct;
SCT(1).Name='DavisStr';
SCT(1).IJ=[462   671
   559   660];


SCT(2).Name='DenmarkStr';
SCT(2).IJ=[781   608
   839   533];

SCT(3).Name='FramStr';
fmt='/home/ddmitry/codes/anls_mtlb_utils/hycom_NOPP_rivers/anls_fluxes_straits/Fram_Sct_mooring008.mat';
IJ=load(fmt); % Fram - close to moorings
SCT(3).IJ=[908  904; ...
           1083 939];
SCT(3).I=IJ.Isct;
SCT(3).J=IJ.Jsct;

SCT(4).Name='NaresStr';
SCT(4).IJ=[756        1062
         756        1083];

SCT(5).Name='BarentsOp';
SCT(5).IJ=[1141         886
        1187         836
        1242         718];

nsct=length(SCT);

% orientation points to determine 
% +/- fluxes across the sections
% Poisitive is northward or towards the AO
IJPR=[636        1209
        1177        1105
         979        1229
        1147        1067
        1296        1055];

ni=length(IJPR);
for ii=1:ni
  i1=IJPR(ii,1);
  j1=IJPR(ii,2);
  IPR(ii)=sub2ind(size(HH),j1,i1);
end

for ii=1:nsct
  SCT(ii).IJPR(:)=IJPR(ii,:);
  SCT(ii).IPR=IPR(ii);
end
%keyboard
% Change orientation of segments if needed
% segments go from i1 to i2 and i1<i2
% in Nares Str: j1 to j2, j1<j2
P=[0, 1; 1, 0]; % permutation matrix
for ip=1:nsct
  i1=SCT(ip).IJ(1,1);
  i2=SCT(ip).IJ(2,1);
  j1=SCT(ip).IJ(1,2);
  j2=SCT(ip).IJ(2,2);

  if i1>i2
    dmm=SCT(ip).IJ;
    dmm=P*dmm;
    SCT(ip).IJ=dmm;
  elseif i1==i2 & j1>j2
    dmm=SCT(ip).IJ;
    dmm=P*dmm;
    SCT(ip).IJ=dmm;
  end
end

for ip=1:nsct
  fprintf('Section: %s\n',SCT(ip).Name);
  IJs=SCT(ip).IJ;

  if ~strncmp(SCT(ip).Name,'Fram',4)
				nij=size(IJs,1);
				IIs=[];
				JJs=[];
				for ii=1:nij-1
						i1=IJs(ii,1);
						i2=IJs(ii+1,1);
						j1=IJs(ii,2);
						j2=IJs(ii+1,2);
						[I,J]=sub_xsct_indx(i1,j1,i2,j2);
						if size(I,1)==1;
								I=I';
								J=J';
						end
			
      if isempty(IIs)	
						  IIs=[IIs;I];
						  JJs=[JJs;J];
      else
        IIs=[IIs;I(2:end)];
        JJs=[JJs;J(2:end)];
      end
				end;

				SCT(ip).I=IIs;
				SCT(ip).J=JJs;
  end
  
  IIs=SCT(ip).I;
  JJs=SCT(ip).J;
  nsg=length(IIs);
  clear XX YY
  for ii=1:nsg
    i0=IIs(ii);
    j0=JJs(ii);
    XX(ii)=LON(j0,i0);
    YY(ii)=LAT(j0,i0);
  end

  SCT(ip).long=XX;
  SCT(ip).latd=YY;
end




return

