function [PRTCL,UN,VN]=sub_seed_particles(HH,LON,LAT,fina,finb,plevel);
% Seed 1 prtcl per 1 grid pnt
% at level = 
ni=3;
iBG = [536, 557+ni;...  % 557+6
       460, 569+ni;...  % 569+6
       460, 569-ni;...
       536, 557-ni];

[mm,nn]=size(LON);
UN=[];
VN=[];
[XX,YY] = meshgrid((1:nn),(1:mm));
INP = inpolygon(XX,YY,iBG(:,1),iBG(:,2));
IN = find(INP==1);
[JJ,II]=ind2sub(size(HH),IN);

fprintf('===> Initial seeding: %i particles\n\n',length(II));
chck=0;
if chck>0
  figure(10); clf;
  contour(HH,[0 0],'k');
  hold on;
  contour(HH,[-1000 -1000],'b');
  contour(HH,[-100 -100],'c');
  plot(iBG(:,1),iBG(:,2),'r.-');
%  contour(LON,[-153 -153],'g');
%  contour(LON,[-130 -130],'g');
%  contour(LAT,[73 73],'g');
%  contour(LAT,[83 83],'g');
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

if ~isempty(fina)
%[F,n,m,l] = read_hycom(fina,finb,'u-vel.');
  [F,n,m,l] = read_hycom(fina,finb,'u-vel.','r_layer',plevel);
  F(F>1e6)=nan;
  UN=squeeze(F);
  [F,n,m,l] = read_hycom(fina,finb,'v-vel.','r_layer',plevel);
  F(F>1e6)=nan;
  VN=squeeze(F);
end

nprt=0; % particle counter
  
PRTCL.TRACK(icnt).I = II;
PRTCL.TRACK(icnt).J = JJ;



return