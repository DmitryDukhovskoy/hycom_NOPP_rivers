function x0 = sub_RegulaFalsi(fn,a0,b0,t0);
% Solve scalar function f(x) = 0
% using Regula Falsi iteration
% Similar to secant method but
% keep two points (x(k)) and x(l) such that
% f(x(k))*f(x(l)) < 0, where k and l are iteration steps
% i.e. always stay on the curve (f(x)) such that 1 point is >0 and 1 is <0
%
% a0, b0 - points of the iteration interval x=[a0,b0], root is in this interval
% t0 - time when GFWA is estimated (yrs, e.g. 24 yrs since 1993)
% fn - function
eps0 = 1e-3;
fprintf('Regula Falsi iteration function %i\n',fn);

if fn==1
  fa = function1(a0,t0);
  fb = function1(b0,t0); 
elseif fn==2
  fa = function2(a0,t0);
  fb = function2(b0,t0); 
elseif fn==3
  fa = function3(a0,t0);
  fb = function3(b0,t0);
end

q1=(fb-fa)/(b0-a0);
xk=a0-fa/q1;
if fn==1
  fxk = function1(xk,t0);
elseif fn==2
  fxk = function2(xk,t0);
elseif fn==3
  fxk = function3(xk,t0);
end

ak = a0;
bk = b0;
fak = fa;

clear X
k=0;
while (abs(fxk))>eps0
  k=k+1;
  if fxk*fak<0 
    ak1=ak;
    bk1=xk;
  else
    ak1=xk;
    bk1=bk;
  end

  if fn==1
    fak1 = function1(ak1,t0);
    fbk1 = function1(bk1,t0);
  elseif fn==2
    fak1 = function2(ak1,t0);
    fbk1 = function2(bk1,t0);
  elseif fn==3
    fak1 = function3(ak1,t0);
    fbk1 = function3(bk1,t0);
  end
  qk = (fbk1-fak1)/(bk1-ak1);
  xk = ak1-fak1/qk;
  if fn==1
    fxk = function1(xk,t0);
  elseif fn==2
    fxk = function2(xk,t0);
  elseif fn==3
    fxk = function3(xk,t0);
  end
  ak=ak1;
  bk=bk1;
  fak=fak1;
  fbk=fbk1;
  X(k)=xk;

  if k>1e6,
    error('Regular Falsi not convergent ...');
  end
end

x0 = X(k);
fprintf('RF iteration converged in %i steps, eps0=%8.6d, root=%10.5f\n\n',k,eps0,x0);

return
