function v=fzerov(f,fv,u0,tol)
%Finding zero of f(u0,v) with given initial guess v=0

%Default absolute tolerance
if nargin==3
    tol=eps;
end
[m,n]=size(u0);
u0=u0(:);
v0=zeros(m*n,1);
v=v0;
errf=abs(f(u0,v0));
errv=ones(m*n,1);
j=1:m*n;
while max(errf)>tol && max(errv)~=0
    v=v0-f(u0,v0)./fv(u0,v0); %Newton iteration
    errf=abs(f(u0,v0));
    errv=abs(v-v0);
    v0=v;
end
v=reshape(v,m,n);