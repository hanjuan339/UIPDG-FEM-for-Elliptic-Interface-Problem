function q = quad1(fun,a,b,p,cv)
%Numerically evaluate integral for p-th order UIPFEM.

%   Q = QUAD1(FUN,A,B,CV,P) approximates the integrals of FUN(X) over intervals
%   A <= X <= B, where A and B may be column vectors of the same size. CV approximates the angles
%   between normal vectors at the end points of interface segement(s).
%   FUN is a function handle.

n=size(a,1);
if n==0
    q=[];
    return;
end

if nargin==4
    cv=zeros(n,1);
end

q=zeros(n,1);
for j=1:2
    if j==1
        l=cv<0.05;
        m=2*p*(p+1)+1;%2*p*(p+1)+1;%2*p;
    else
        l=~l;
        m=2*p*(p+1)+1;%2*p*(p+1)+1;%4*p;
    end
    if nnz(l)
        [x,w]=lgwt(m,0,1);%Gaussian nodes and weights in [0,1]
        X=repmat(a(l),1,m)+(b(l)-a(l))*x';
        F=zeros(size(X));
        for j=1:m
            F(:,j)=fun(X(:,j),l);
        end
        q(l)=(b(l)-a(l)).*(F*w);
    end
end
