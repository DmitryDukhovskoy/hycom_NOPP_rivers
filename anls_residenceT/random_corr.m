n=1000;
X=rand(n,1);
Y=rand(n,1);
X=X-mean(X);
Y=Y-mean(Y);
Rho=xcorr(X,Y,'coeff',3);
rho = mean(X.*Y)./(std(X).*std(Y))

t=[1:n]';
alf=0.0015;
bet=0.0012;
X1=X+alf.*t;
Y1=Y+bet*t;
X1=X1-mean(X1);
Y1=Y1-mean(Y1);
rho1 = mean(X1.*Y1)./(std(X1).*std(Y1))


figure(1); clf;
subplot(2,1,1);
hold on;
plot(X);
plot(Y);
stl=sprintf('Corr=%5.3f',rho);
title(stl);
set(gca,'tickdir','out',...
        'fontsize',14);

subplot(2,1,2);
hold on;
plot(X1);
plot(Y1);
stl=sprintf('Corr=%5.3f',rho1);
title(stl);
set(gca,'tickdir','out',...
        'fontsize',14);




