function phi=upfemrefbas(x,y,i,j,p,xn)
%The (i,j) basis function on refernce element
%   p is the degree of element polynomials (in each x and y, respectively)
%   xn: nodal points in [0,1].      

[m,n]=size(x);
x=reshape(x,1,m*n);
y=reshape(y,1,m*n);

xxi=repmat(x,p+1,1)-repmat(xn,1,m*n);
xxi(i,:)=[];
numx=prod(xxi,1); %Numerator of the part of the basis function related to x
denx=xn(i)-xn;
denx(i)=[];
denx=prod(denx); %Denominator of the part of the basis function related to x

yyj=repmat(y,p+1,1)-repmat(xn,1,m*n);
yyj(j,:)=[];
numy=prod(yyj,1); %Numerator of the part of the basis function related to y
deny=xn(j)-xn;
deny(j)=[];
deny=prod(deny); %Denominator of the part of the basis function related to y

phi=(numx/denx).*(numy/deny);
phi=reshape(phi,m,n);
