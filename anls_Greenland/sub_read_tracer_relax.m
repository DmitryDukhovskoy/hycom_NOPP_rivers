function F = sub_read_tracer_relax(ftrca,ftrcb,IND,imonth,nlev);
% Read tracer concentration - relax file
% used for the simulation
% this file is created during tracer or river preparation

%regn='ARCc0.08';
%ntopo=9;      % HYCOM topo version
kplt=1; % layer to plot

inc1=IND.i1;
inc2=IND.i2;
jnc1=IND.j1;
jnc2=IND.j2;
djnc=IND.dj;
dinc=IND.di;


fidb = fopen(ftrcb,'r');  % read I,J from *.b
for nl=1:5
  aa=fgetl(fidb);
  disp(aa);
end

is=strfind(aa,'= ');
[ID,JD] = strread(aa(is+1:end),'%d%d');

disp(['Grid I=',num2str(ID),' J=',num2str(JD)]);
IJDM=ID*JD;
npad=4096-mod(IJDM,4096);
toto=ones(npad,1);
fida = fopen(ftrca,'r');  % 
for im=12:12
  fprintf('Tracer input: Reading mo=%i lev=%i\n',im,kplt);
  nrec=(im-1)*nlev+kplt-1;
  stat=fseek(fida,nrec*(IJDM+npad)*4,-1);
  dmm=fread(fida,IJDM,'float32','ieee-be');  % read 2D field (1 layer)
%  dm1=fread(fid1,npad,'float32','ieee-be');  % read npad 

  A=reshape(dmm,ID,JD);
  A=A';
  
  ll=find(A>0);
  fprintf('# of river cells %i\n',length(ll));

  mnA=min(min(A(A>0)));
  mxA=max(max(A));
  
  F=A(jnc1:jnc2,inc1:inc2);
  F(F==0)=nan;
  
  szz = sprintf('min C=%5.2d, \nmax C=%5.2d',mnA,mxA);
  fprintf('%s\n',szz);
  
%  clf;
%  contour(H,[0 0],'Color',[0.8 0.8 0.8]);
%  hold on;
%  pcolor(F); shading flat;
%  caxis([0 100]);
%  axis('equal');
%  set(gca,'xlim',[250 450],'ylim',[250 550]);
%  colorbar;
%  spp=sprintf('Tracer conc., mo=%i, level=%i', im,kplt);
%  title(spp);
%  text(255,500,szz);
%  drawnow

end;




return