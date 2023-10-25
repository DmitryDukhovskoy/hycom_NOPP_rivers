% Plot time series of fluxes
%
function sub_plotflx(Vflx, nfg, sunits, tk1, dtk, tk2, TM, clr1, clr2, flxnm);

DV = datevec(TM);
% Monthly climatology
% Group by months
for im=1:12
  I=find(DV(:,2)==im);

  dmm=Vflx(I);
  F.Vmn(im)=nanmean(dmm);
  F.Vlw(im)=prctile(dmm,25);
  F.Vup(im)=prctile(dmm,75);
end

figure(nfg); clf;
set(gcf,'Position',[1378 532 1171 784]);
axes('Position',[0.1 0.55 0.83 0.38]);
%axes('Position',POS(2,:));
hold on

vflx=F.Vmn; 
v1=F.Vlw;
v2=F.Vup;

if isempty(clr1), clr1=[0 0.4 0.8]; end
if isempty(clr2), clr2=[0 0.6 1]; end;

mnth=[1:12];
plot(mnth,vflx,'Linewidth',2.5,'Color',clr1); % mSv, + to the Arctic

for ik=1:12
  plot([ik ik],[v1(ik) v2(ik)],'-','Color',clr2);
  plot([ik-0.1 ik+0.1],[v1(ik) v1(ik)],'-','Color',clr2);
  plot([ik-0.1 ik+0.1],[v2(ik) v2(ik)],'-','Color',clr2);
end
% Overall statistics:
mVflx=mean(Vflx);
vlw=prctile(Vflx,25);
vup=prctile(Vflx,75);

stxt=sprintf('%s: %5.2f IQR: %5.2f/%5.2f %s',flxnm,mVflx,vlw,vup,sunits);
title(stxt);

dyy=0.03*abs(max(v2)-min(v1));
yl1=min(v1)-dyy;
yl2=max(v2)+dyy;
set(gca,'tickdir','out',...
 'xlim',[0.5 12.5],...
 'xtick',[1:12],...
 'ylim',[yl1 yl2],...
 'ytick',[tk1:dtk:tk2],...
 'xgrid','on',...
 'ygrid','on',...
 'Fontsize',14);



%
% Monthly means time series:
yr1=DV(1,1);
yr2=DV(end,1);
icc=0;
jjr=0;
for iyr=yr1:yr2
  for im=1:12
    icc=icc+1;
    I=find(DV(:,1)==iyr & DV(:,2)==im);

    dmm=Vflx(I);
    F.Vmo(icc,1)=mean(dmm);
    F.VMlw(icc,1)=prctile(dmm,25);
    F.VMup(icc,1)=prctile(dmm,75);

    tmo(icc,1)=iyr+(im-1)/12;
  end

% Annual mean:
  jjr=jjr+1;
  J=find(DV(:,1)==iyr);
  dmm=Vflx(J);
  F.Vyr(jjr)=mean(dmm);
  
end


%yl1=floor(prctile(Vflx,1));
%yl2=ceil(prctile(Vflx,99));
dyy=0.03*abs(max(F.VMup)-min(F.VMlw));
yl1=min(F.VMlw)-dyy;
yl2=max(F.VMup)+dyy;
axes('Position',[0.1 0.09 0.83 0.38]);
%axes('Position',POS(2,:));
hold on
plot(tmo,F.Vmo,'Linewidth',2.5,'Color',clr1);
plot(tmo,F.VMlw,'Linewidth',1.5,'Color',clr2);
plot(tmo,F.VMup,'Linewidth',1.5,'Color',clr2);

for iyr=yr1:yr2
  plot([iyr iyr],[yl1 yl2],'--','Color',[0.8 0.8 0.8]);
end

% Plot annual mean:
jjr=0;
for iyr=yr1:yr2
  jjr=jjr+1;
  vyr=F.Vyr(jjr);
  plot([iyr iyr+1],[vyr vyr],'-','Linewidth',1.6,'Color',[0.5 0.5 0.5]);
end

%keyboard
set(gca,'tickdir','out',...
 'xlim',[yr1 yr2+1],...
 'xtick',[yr1:2:yr2+1],...
 'ylim',[yl1 yl2],...
 'ytick',[tk1:2*dtk:tk2],...
 'xgrid','on',...
 'ygrid','on',...
 'Fontsize',14);



return
