function rskp=sub_skip_rivers(rnm);
%
% Skip tributaries and unneeded rivers in UK
rskp=1;

if strncmp(rnm,'Peza',4)
  return;
end
if strncmp(rnm,'Sula',4)
  return;
end
if strcmp(rnm,'Poluy'); % Ob tributary
  return;
end
if strcmp(rnm,'Bolshoy Anyuy'); % Kolyma tributary
  return;
end
if strcmp(rnm,'Tamar'); % 
  return;
end
if strcmp(rnm,'Exe'); % 
  return;
end
if strcmp(rnm,'Taff'); % 
  return;
end
if strcmp(rnm,'Tywi'); % 
  return;
end
if strcmp(rnm,'Conwy'); % 
  return;
end
if strcmp(rnm,'Dee (Afon Dyfrd'); % 
  return;
end
if strcmp(rnm,'Severn'); % 
  return;
end
if strcmp(rnm,'Usk'); % 
  return;
end
if strcmp(rnm,'Wye'); % 
  return;
end
if strcmp(rnm,'Tay'); % 
  return;
end
if strcmp(rnm,'Dee (Royal Dee)'); % 
  return;
end
if strcmp(rnm,'Spey'); % 
  return;
end
if strcmp(rnm,'Tyne'); % 
  return;
end
if strcmp(rnm,'Wharfe'); % 
  return;
end
if strcmp(rnm,'Derwent-EU'); % 
  return;
end
if strcmp(rnm,'Trent'); % 
  return;
end
if strcmp(rnm,'Bedford'); % 
  return;
end
if strcmp(rnm,'Medway'); % 
  return;
end
if strcmp(rnm,'Great Stour'); % 
  return;
end
if strcmp(rnm,'Itchen'); % 
  return;
end
if strcmp(rnm,'Kennet'); % 
  return;
end
if strcmp(rnm,'Bedford Ouse'); % 
  return;
end
if strcmp(rnm,'Avon (Bristol A'); % 
  return;
end

rskp = 0;

return