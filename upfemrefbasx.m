function phix=upfemrefbasx(x,y,i,j,p,xn)
%The partial derivative of (i,j) basis function on refernce element w.r.t x
%   p is the degree of element polynomials (in each x and y, respectively).
%   1<=i,j<=p+1.
%   xn: nodal points in [0,1].      

[m,n]=size(x);
x=reshape(x,1,m*n);
y=reshape(y,1,m*n);

xxi=repmat(x,p+1,1)-repmat(xn,1,m*n);
xxi(i,:)=[];
denx=xn(i)-xn;
denx(i)=[];
denx=prod(denx); %Denominator of the part of the basis function related to x

yyj=repmat(y,p+1,1)-repmat(xn,1,m*n);
yyj(j,:)=[];
numy=prod(yyj,1); %Numerator of the part of the basis function related to y
deny=xn(j)-xn;
deny(j)=[];
deny=prod(deny); %Denominator of the part of the basis function related to y

%dphi/dx
if p==1
    numdx=1;
else
    numdx=0;
    for k=1:p
        xx=xxi;
        xx(k,:)=[];
        numdx=numdx+prod(xx,1);
    end
end
phix=(numdx/denx).*(numy/deny);
phix=reshape(phix,m,n);
