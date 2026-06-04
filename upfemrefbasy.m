function phiy=upfemrefbasy(x,y,i,j,p,xn)
%The partial derivative of (i,j) basis function on refernce element w.r.t y
%   p is the degree of element polynomials (in each x and y, respectively).
%   1<=i,j<=p+1.
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
deny=xn(j)-xn;
deny(j)=[];
deny=prod(deny); %Denominator of the part of the basis function related to y

%dphi/dy
if p==1
    numdy=1;
else
    numdy=0;
    for k=1:p
        yy=yyj;
        yy(k,:)=[];
        numdy=numdy+prod(yy,1);
    end
end
phiy=(numx/denx).*(numdy/deny);
phiy=reshape(phiy,m,n);
