function Rnew = sub_rivers2coast(HH,R);
% Adjust river sources (R) to the new coast line
% in the new Topo (HH)
% old topo = HHo

Rtot1=nansum(nansum(R));

% Onland river sources:
% that need to be moved
IR=find(R>0 & HH>=0);
nr=length(IR);
nr0=nr; 

Inl=find(R>0 & HH<0);
%Rnew = HH*0;
%Rnew = R(Inl);
Rnew = R;

fprintf(':: Found %i river grid cells on land\n',nr);
if nr==0; return; end;

[mm,nn]=size(R);

% identify a cluster of points
% around a river source
kk=1;
cc=0;
di=20; % search radius for river points
dj=di;
while nr>0
  ia=min(find(~isnan(IR)));
  if isempty(ia); break; end;
  cc=cc+1;
  if cc>nr0,
    fprintf('sub_rivers2coast: endless loop');
    error('check river points cc=%i\n',cc);
  end

  % Define a subset with the river source
  % and work with the subset
  ir0=IR(ia);
  [j0,i0]=ind2sub(size(R),ir0);
% find closest river points
  j1=j0-dj;
  ji=max([1,j1]);
  j2=j0+dj;
  j2=min([mm,j2]);
  i1=i0-di;
  i1=max([1,i1]);
  i2=i0+di;
  i2=min([nn,i2]);
% subsets for given river:  
  sr=R(j1:j2,i1:i2);
  sh=HH(j1:j2,i1:i2);
  isb=find(sr>0);
  RG=sr(isb);
  RG=sort(RG,'descend');
% Start from max runoff location
% move it to the ocean if needed
  [jmx,imx]=find(sr==RG(1)); % location of the max runoff
  hmx=sh(jmx,imx);
  jm0=jmx;
  im0=imx;
  if hmx>=0,
    [Joc,Ioc]=find(sh<0);
    do=sqrt((Joc-jmx).^2+(Ioc-imx).^2);
    idd=find(do==min(do),1);
    jm0=Joc(idd);
    im0=Ioc(idd);
  end

% Distribute runoff in the ocean grid cells
% closest to the max runoff point
  srn=sh*0;
  srn(sh>=0)=9e10; % fill land
  for kr=1:length(RG);
    srn(jm0,im0)=RG(kr);
    [je,ie]=find(srn==0);
%    dd=sqrt((je-jm0).^2+(ie-im0).^2);
    dx=sqrt((je-jmx).^2+(ie-imx).^2); % also closest to the old max
%    idd=find(dd==min(dd));
%    if length(idd)>1
%      ip=find(dx(idd)==min(dx(idd)));
%      idd=idd(ip);
%    end
    idd=find(dx==min(dx),1);   
    jm0=je(idd);
    im0=ie(idd);
  end
  srn(srn>1e10)=nan;
  srn(srn==0)=nan;
    
  
  
  rtot1=nansum(nansum(sr));
  rtot2=nansum(nansum(srn));
  
  if abs(1-rtot1/rtot2)>0.001
    fprintf('Check total runoff in relocated river\n');
    fprintf('Old %8g, New %8g\n',rtot1,rtot2);
    keyboard;
  end
  
  Rnew(j1:j2,i1:i2)=srn;
%  R(j1:j2,i1:i2)=0;
  
  IR=find(Rnew>0 & HH>=0);
  nr=length(IR);
  
end;

Rnew(isnan(Rnew))=0; % NO LAND mask in rivers

Rtot2=nansum(nansum(Rnew));
  if abs(1-Rtot1/Rtot2)>0.001
    fprintf('Mismatch: domain-integrated River runoff in corrected data\n');
    fprintf('Old %8g, New %8g\n',Rtot1,Rtot2);
    keyboard;
  end

fprintf(':: Rivers fixed\n');

return