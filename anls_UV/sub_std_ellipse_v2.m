function sub_std_ellipse_v2(uu,vv,ii,jj,sclv,v_col,cf,lwd,beta);
% draw std ellipse
%
um2cm=1;
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
ss=sqrt(mux*mux+muy*muy);


%if isempty(sclv), sclv=20; end;
%sclv=20;   % scaling=# of gr.pnts/1cm/s of the velocity, ellipse
%cf=0.35;
%beta=20;
col=v_col;
%lwd=2.;
[Xh,Yh]=get_arrowF(ii,ii+mux*sclv,jj,jj+muy*sclv,cf,beta,col,lwd);
x0=Xh(1);
y0=Yh(1);

%keyboard
%      plot([ii ii+mux],[jj jj+muy],'r');
pp8=plot(ii,jj,'k.','MarkerSize',9);
%set(pp8,'Color',col);
%      hold on
col_ell=[0 0. 0.4];   % ellipse color
lwde=1.5;
plot([x0-mj(1)*sclv x0+mj(1)*sclv],...
     [y0-mj(2)*sclv y0+mj(2)*sclv],...
     'Color',col_ell,'LineWidth',lwde);
plot([x0-mn(1)*sclv x0+mn(1)*sclv],...
     [y0-mn(2)*sclv y0+mn(2)*sclv],...
     'Color',col_ell,'LineWidth',lwde);
     
% Plot ellipse:
tt=(0:0.1:2*pi+0.1);
aa=sqrt((mj*sclv)*(mj*sclv)');
bb=sqrt((mn*sclv)*(mn*sclv)');
A=[aa*cos(tt);bb*sin(tt)];
alf=atan2(mj(2),mj(1));
Rot=[cos(-alf), sin(-alf); -sin(-alf), cos(-alf)];
A2=Rot*A;
Ell=[A2(1,:)+x0;A2(2,:)+y0];
plot(Ell(1,:),Ell(2,:),'Color',col_ell,'LineWidth',lwde)





return