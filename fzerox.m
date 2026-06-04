function x=fzerox(f,fx,x0,y0,tol)
%Finding zero of f(x,y0) with given initial guess x0

%Default absolute tolerance
if nargin==4
    tol=eps;
end
[m,n]=size(y0);
y0=y0(:);
x0=x0(:);
x0=x0+rand(size(x0))*eps;
if length(x0)==1
    x0=repmat(x0,m*n,1);
end
x=x0;
errf=abs(f(x0,y0));
errx=ones(m*n,1);
j=1:m*n;
while max(errf(j))>tol & max(errx(j))~=0
    x(j)=x0(j)-f(x0(j),y0(j))./fx(x0(j),y0(j)); %Newton iteration
    errf(j)=abs(f(x(j),y0(j)));
    errx(j)=abs(x(j)-x0(j));
    x0(j)=x(j);
    j=errf>tol & errx~=0;
end
x=reshape(x,m,n);