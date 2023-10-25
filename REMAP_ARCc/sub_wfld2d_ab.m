function sub_wfld2d_ab(nlev,F,IDM,JDM,fida_NEW,fidb_NEW,cline,TDENS);
% Continous writing of 2D fields
% write 2D fields into *a file and update *b
% modifying string cline, which has a general
% pattern of *b file line
% Field F is [k,JDM,IDM]
% output files already open for writing
% nlev - layer # where to put fields
%
hg=2^100;

F=squeeze(F);
a1=size(F);
f2d=(length(a1)==3); % 2d or 3d array

IJDM=IDM*JDM;
npad=4096-mod(IJDM,4096);
toto=ones(npad,1);

%keyboard
if f2d
  A = squeeze(F(nlev,:,:));
else
  A=squeeze(F);
end

A = A';
dmm = reshape(A,IJDM,1);
I=find(isnan(dmm));
umm=nanmean(dmm);
dmm(I)=hg;

fwrite(fida_NEW,dmm,'float32','ieee-be');
fwrite(fida_NEW,toto,'float32','ieee-be');

% Update *b file:      
Inh=find(dmm<hg/10);
minv=min(dmm(Inh));
maxv=max(dmm(Inh));

ir = strfind(cline,'=');
SD = sscanf(cline(ir+1:end),'%f');
td0=TDENS(nlev);

% Parse string from old nest *b:
aa2 = sub_parse_string_B(cline,nlev,td0,minv,maxv);

% Write *b file:
fprintf('<= %s\n',cline);
fprintf('=> %s\n',aa2);
fprintf(fidb_NEW,[aa2,'\n']);



return