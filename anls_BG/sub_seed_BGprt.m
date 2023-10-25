function [PRTCL,UN,VN]=sub_seed_particles(HH,LON,LAT,fina,finb);
% Seed 1 prtcl per 1 grid pnt
% in the upper layer
iBG = [483, 1532;...
       649, 1380; ...
       732, 1477;...
       601, 1637;...
       483, 1532];

[mm,nn]=size(LON);

[XX,YY] = meshgrid((1:nn),(1:mm));
INP = inpolygon(XX,YY,iBG(:,1),iBG(:,2));
IN = find(INP==1);
[JJ,II]=ind2sub(size(HH),IN(1:5:end));

chck=0;
if chck>0
  figure(10); clf;
  contour(HH,[0 0],'k');
  hold on;
  contour(HH,[-1000 -1000],'b');
  plot(iBG(:,1),iBG(:,2),'r.-');
  contour(LON,[-153 -153],'g');
  contour(LON,[-130 -130],'g');
  contour(LAT,[73 73],'g');
  contour(LAT,[83 83],'g');
  [JJ,II]=ind2sub(size(HH),IN);
  plot(II,JJ,'y.');
  axis('equal');
end

rg=9806;  % convert pressure to depth, m
PRTCL.TRACK.I=[];
PRTCL.TRACK.J=[];
%PRTCL.TRACK.T=[];
%PRTCL.TRACK.S=[];
%PRTCL.TRACK.Z=[];

PRTVRT=[]; % array with particle distribution

icnt=1;
nI=length(IN);

%[F,n,m,l] = read_hycom(fina,finb,'u-vel.');
[F,n,m,l] = read_hycom(fina,finb,'u-vel.','r_layer',1);
F(F>1e6)=nan;
UN=squeeze(F);
[F,n,m,l] = read_hycom(fina,finb,'v-vel.','r_layer',1);
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
  
PRTCL.TRACK(icnt).I = II;
PRTCL.TRACK(icnt).J = JJ;



return