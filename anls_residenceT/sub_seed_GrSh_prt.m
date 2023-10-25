function [PRTCL,UN,VN]=sub_seed_GrSh_prt(slr,HH,LON,LAT,fina,finb,Np0,nsim);
% Seed 1 prtcl per 1 grid pnt
% in the layer slr
iBG = [513   624
   561   618
   567   524
   603   432
   590   401
   532   464
   510   570];

[mm,nn]=size(LON);

[XX,YY] = meshgrid((1:nn),(1:mm));
INP = inpolygon(XX,YY,iBG(:,1),iBG(:,2));
IN = find(INP==1 & (HH>-600 & HH<-50));

% Subset particles to get no more than Np0
np=length(IN);
for ii=1:nsim-1  % dummy random to have different random for nsim
  dmm=rand(Np0,1);
end

if np<=Np0  % 
   fprintf('Requested # of particles %i > grid points: %i\n',Np0,np);
  [JJ,II]=ind2sub(size(HH),IN);
else
  rng=('shuffle');
  dmm=1+round(rand(Np0,1)*(np-1));
  Ipp=unique(dmm);
  cc=0;
% No repetition in indices
  while length(Ipp)<Np0
    b=1+round(rand(1)*(np-1));
    Ipp=[Ipp;b];
    Ipp=unique(Ipp);
    cc=cc+1;
    if cc>1000
      error('sub_seed_GrSh_prt: Endless loop random number generator');
    end
  end
  
  [JJ,II]=ind2sub(size(HH),IN(Ipp));
end

PRTCL.TRACK.I=[];
PRTCL.TRACK.J=[];
%PRTCL.TRACK.T=[];
%PRTCL.TRACK.S=[];
%PRTCL.TRACK.Z=[];
UN =[];
VN = [];

icnt=1;
nI=length(IN);
PRTCL.TRACK(icnt).I = II;
PRTCL.TRACK(icnt).J = JJ;
    
if isempty(fina) | isempty(finb)
  fprintf('Not requested U and V for particles ...\n');
  return
end  
  

%keyboard

chck=0;
if chck>0
  figure(10); clf;
  contour(HH,[0 0],'k');
  hold on;
  contour(HH,[-8000:1000:-1000],'b');
  plot(iBG(:,1),iBG(:,2),'r.-');
  [JJ,II]=ind2sub(size(HH),IN);
  plot(II,JJ,'y.');
  axis('equal');
  set(gca,'xlim',[820 1200],...
	  'ylim',[450 950]);
end

rg=9806;  % convert pressure to depth, m

%PRTVRT=[]; % array with particle distribution


%[F,n,m,l] = read_hycom(fina,finb,'u-vel.');
[F,n,m,l] = read_hycom(fina,finb,'u-vel.','r_layer',slr);
F(F>1e6)=nan;
UN=squeeze(F);
[F,n,m,l] = read_hycom(fina,finb,'v-vel.','r_layer',slr);
F(F>1e6)=nan;
VN=squeeze(F);

%[F,n,m,l] = read_hycom(fina,finb,'thknss');
%F=F./rg;
%F(F>1e10)=nan;
%dH=squeeze(F); 

%[F,n,m,l] = read_hycom(fina,finb,'temp');
%F(F>1e10)=nan;
%T=squeeze(F); 

%[F,n,m,l] = read_hycom(fina,finb,'salin');
%F(F>1e10)=nan;
%S=squeeze(F); 

nprt=0; % particle counter
  



return
