function [FM,TMm] = sub_day2month(F,TM);
%
% From daily fields
% calculate monthly mean fields
%
fprintf('Calculating monthly means ...\n');
%keyboard

DV=datevec(TM);
cp=0;
clear FM
for iyr=DV(1,1):DV(end,1),
  IY = find(DV(:,1)==iyr);
  dv = DV(IY,:);
  ff = F(IY,:);
  for im=1:12
    IM = find(dv(:,2)==im);
    cp=cp+1;
    FM(cp,:) = mean(ff(IM,:),1);
    TMm(cp,1) = datenum(iyr,im,15);
  end
end

return