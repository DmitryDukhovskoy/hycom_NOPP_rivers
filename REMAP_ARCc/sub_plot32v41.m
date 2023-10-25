function sub_plot32v41(Fnew,Fold,HH,IJ,ZZn,ZZo,vname);
% Check interpolated fields
% 32 -> 41 v. layers
% F - interpolated field
% HH - bathymetry, should be same in old(32) and new(41)
%     fondigurations
% IJ - indices of x-section
% ZZn, ZZo - new and old isopycnal interface depths
%

jj1=IJ(1,2);
jj2=IJ(2,2);
ii1=IJ(1,1);
ii2=IJ(2,1);

A=squeeze(Fnew(:,jj1:jj2,ii1:ii2));
[lr41,np]=size(A);
% Add extra layer:
A(lr41+1,:)=A(lr41,:);


Zn=squeeze(ZZn(:,jj1:jj2,ii1:ii2));
for kk=1:lr41+1
  XLn(kk,:)=[1:np];
end

figure(10); clf;
axes('Position',[0.08 0.52 0.9 0.42]);
pcolor(XLn,Zn,A); shading flat;
hold on;
set(gca,'Color',[0 0 0]);
colorbar
stt=sprintf('Interpolated %s, I: %i-%i, J: %i-%i',...
	    vname,ii1,ii2,jj1,jj2);
title(stt);


axes('Position',[0.08 0.05 0.9 0.42]);
pcolor(XLn,Zn,A); shading flat;
hold on;
set(gca,'Color',[0 0 0]);
colorbar

lr00=24;
for kk=1:lr41
  zz=Zn(kk,:);
  xx=XLn(kk,:);
  if mod(kk,5)==0;
    plot(xx,zz,'-','Color',[0.4 0.4 0.4]);
  else
    plot(xx,zz,'b:','Color',[0.6 0.6 0.6]);
  end 
  if kk==lr00+1
    plot(xx,zz,'g--');
  end
  
end

% Plot old field:
Ao=squeeze(Fold(:,jj1:jj2,ii1:ii2));
[lr32,np]=size(Ao);
% Add extra layer:
Ao(lr32+1,:)=Ao(lr32,:);

Zo=squeeze(ZZo(:,jj1:jj2,ii1:ii2));
for kk=1:lr32+1
  XLo(kk,:)=[1:np];
end
%txtb='ddmitry@mars: hycom_NOPP_rivers/REMAP_ARCc/remap32to41lrs_fixedZ.m';
txtb='ddmitry@mars: hycom_NOPP_rivers/REMAP_ARCc/sub_plot32v41.m';
bottom_text(txtb);



figure(11); clf;
axes('Position',[0.08 0.52 0.9 0.42]);
pcolor(XLo,Zo,Ao); shading flat;
hold on;
colorbar
set(gca,'Color',[0 0 0]);
stt=sprintf('Original %s, I: %i-%i, J: %i-%i',...
	    vname,ii1,ii2,jj1,jj2);
title(stt);


axes('Position',[0.08 0.05 0.9 0.42]);
pcolor(XLo,Zo,Ao); shading flat;
hold on;
set(gca,'Color',[0 0 0]);
colorbar

lr00=14;
for kk=1:lr32
  zz=Zo(kk,:);
  xx=XLo(kk,:);
  if mod(kk,5)==0;
    plot(xx,zz,'-','Color',[0.4 0.4 0.4]);
  else
    plot(xx,zz,'b:','Color',[0.6 0.6 0.6]);
  end 
  if kk==lr00+1
    plot(xx,zz,'g--');
  end
  
end
txtb='ddmitry@mars: hycom_NOPP_rivers/REMAP_ARCc/sub_plot32v41.m';
bottom_text(txtb);






return