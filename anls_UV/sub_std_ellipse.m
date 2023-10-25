function sub_std_ellipse(uu,vv,ii,jj,um2cm,sclv);
% draw std ellipse
%
if isempty(um2cm), um2cm=100; end; % m/s-> cm/s
uu=uu*um2cm;
vv=vv*um2cm;
UV(:,1) = uu-mean(uu);
UV(:,2) = vv-mean(vv);
Cv=1/length(UV)*UV'*UV;

% Find eigenvectors (Ev) and eigenvalues (El):
[Ev El]=eig(Cv);
e1=Ev(:,1)';
e2=Ev(:,2)';
lmb1=El(1,1);
lmb2=El(2,2);

% Major axis (Applied multivariate statistical analysis, p. 463):
mj=(lmb1)^0.5*e1;
mn=(lmb2)^0.5*e2;

mux=mean(uu);
muy=mean(vv);

if isempty(sclv), sclv=20; end;
%sclv=20;   % scaling=# of gr.pnts/1cm/s of the velocity, ellipse
cf=0.35;
beta=20;
col=[0 0 0];
lwd=2.;
draw_arrowF(ii,ii+mux*sclv,jj,jj+muy*sclv,cf,beta,col,lwd);
%keyboard
%      plot([ii ii+mux],[jj jj+muy],'r');
pp8=plot(ii,jj,'k.','MarkerSize',9);
%set(pp8,'Color',col);
%      hold on
col_ell=[1 0.3 0];   % ellipse color
lwde=1.5;
plot([ii+(mux-mj(1))*sclv ii+(mux+mj(1))*sclv],...
     [jj+(muy-mj(2))*sclv jj+(muy+mj(2))*sclv],...
     'Color',col_ell,'LineWidth',lwde);
plot([ii+(mux-mn(1))*sclv ii+(mux+mn(1))*sclv],...
     [jj+(muy-mn(2))*sclv jj+(muy+mn(2))*sclv],...
     'Color',col_ell,'LineWidth',lwde);
% Plot ellipse:
tt=(0:0.1:2*pi+0.1);
aa=sqrt((mj*sclv)*(mj*sclv)');
bb=sqrt((mn*sclv)*(mn*sclv)');
A=[aa*cos(tt);bb*sin(tt)];
alf=atan2(mj(2),mj(1));
Rot=[cos(-alf), sin(-alf); -sin(-alf), cos(-alf)];
A2=Rot*A;
Ell=[A2(1,:)+ii+mux*sclv;A2(2,:)+jj+muy*sclv];
plot(Ell(1,:),Ell(2,:),'Color',col_ell,'LineWidth',lwde)





return