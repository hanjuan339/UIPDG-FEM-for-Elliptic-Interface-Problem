function fxy=upfemeval(x,y,lvsfun,fun1,fun2)
%Evaluating piecewise function.
%   fxy=f(x,y) where f=fun1 in Omega_1 and f=fun2 in Omega2.
[m1,m2]=size(x);
x=x(:);
y=y(:);
phi=lvsfun(x,y);
i1=find(phi<=0);
i2=find(phi>0);
fxy=zeros(m1*m2,1);
fxy(i1)=fun1(x(i1),y(i1));
fxy(i2)=fun2(x(i2),y(i2));
fxy=reshape(fxy,m1,m2);
