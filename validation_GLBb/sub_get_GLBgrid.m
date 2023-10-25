function sub_get_GLBgrid(pthglb,expt,fmat);
% Subset topo/grid for
% specified region
% from GLBb grid & topo
fprintf('Subsetting region from GLBb ...\n');
% Indices for N.Atl. subregion:
IND.i1=2774;
IND.i2=3697;
IND.j1=2204;
IND.j2=3197;
i1=IND.i1;
i2=IND.i2;
j1=IND.j1;
j2=IND.j2;

if expt==190 | expt==191
  TV='07';
end


frga   = sprintf('%sregional.grid.a',pthglb);
frgb   = sprintf('%sregional.grid.b',pthglb);
fdptha = sprintf('%sdepth_GLBb0.08_%s.a',pthglb,TV);
fdpthb = sprintf('%sdepth_GLBb0.08_%s.b',pthglb,TV);

fidRGb = fopen(frgb,'r');  % read I,J from regional.grid.b
aa  = fgetl(fidRGb);
dmm = aa(2:8);
IDM = str2num(dmm);
aa = fgetl(fidRGb);
dmm = aa(2:8);
JDM = str2num(dmm);
IJDM = IDM*JDM;
fclose(fidRGb);
npad=4096-mod(IJDM,4096);

fprintf('Topo: IDM=%i, JDM=%i\n',IDM,JDM);

% read lon/lat from GLBb regional grid file
fidRGa = fopen(frga,'r');
[plon,count] = fread(fidRGa,IJDM,'float32','ieee-be');
fseek(fidRGa,4*(npad+IJDM),-1);
[plat,count] = fread(fidRGa,IJDM,'float32','ieee-be');

fprintf('Reading lat/lon for GLBb0.08 - %3.3i ...\n',expt)
plon=(reshape(plon,IDM,JDM))';
plat=(reshape(plat,IDM,JDM))';
fclose(fidRGa);

fdptha = sprintf('%sdepth_GLBb0.08_07.a',pthglb);
fdpthb = sprintf('%sdepth_GLBb0.08_07.b',pthglb);
fida = fopen(fdptha);
HHg = fread(fida,IJDM,'float32','ieee-be');
HHg = (reshape(HHg,IDM,JDM))';
fclose(fida);

I=find(HHg>1e10);
HHg=-HHg;
HHg(I)=100;

LON = plon(j1:j2,i1:i2);
LAT = plat(j1:j2,i1:i2);
HH  = HHg(j1:j2,i1:i2);

Lmx=max(max(LON));
cc=0;
while(Lmx>180);
  cc=cc+1;
  if (cc>100), error('sub_get_GKBgrid: endless loop, LON>180'); end;
  I=find(LON>180);
  LON(I)=LON(I)-360;
  Lmx=max(max(LON));
end;
  
Lmn=min(min(LON));
cc=0;
while(Lmn<-180);
  cc=cc+1;
  if (cc>100), error('sub_get_GKBgrid: endless loop, LON<-180'); end;
  I=find(LON>180);
  LON(I)=LON(I)+360;
  Lmn=min(min(LON));
end;
  

%fmat = sprintf('%sGLBb_T%s_subset_Natl',pthmat,TV);
fprintf('Saving subregion NAtl., T=%s, lon/lat/H\n',TV);
fprintf('%s\n',fmat);
save(fmat,'LON','LAT','HH','IND');

return