function y=fzeroy(f,fy,x0,y0,tol)
%Finding zero of f(x0,y) with given initial guess y0

%Default absolute tolerance
if nargin==4
    tol=eps;
end
[m,n]=size(x0);
x0=x0(:);
y0=y0(:);
y0=y0+rand(size(y0))*eps;
if length(y0)==1
    y0=repmat(y0,m*n,1);
end
y=y0;
errf=abs(f(x0,y0));
erry=ones(m*n,1);
j=1:m*n;
while max(errf(j))>tol & max(erry(j))~=0
    y(j)=y0(j)-f(x0(j),y0(j))./fy(x0(j),y0(j)); %Newton iteration
    errf(j)=abs(f(x0(j),y(j)));
    erry(j)=abs(y(j)-y0(j));
    y0(j)=y(j);
    j=errf>tol & erry~=0;
end
y=reshape(y,m,n);