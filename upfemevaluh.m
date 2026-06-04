function U=upfemevaluh(x,y,u,lvsfun,el1,el2,eli,ltgmap1,ltgmap2,ltgmapi1,ltgmapi2,n,p,xn)
%Evaluating the UPFE soltion 
%   x,y: the coordinates of the points (in the unit square). x,y may be
%   matrics of the same size.
%   u: Nodal values of the discrete solution.
%   lvsfun: the level set function.
%   xn: nodal points in [0,1].      

h=1/n;
ldofs=(p+1)^2;
[row,col]=size(x);
X=x(:);
Y=y(:);
I=max(0,floor(n*(X-eps))); %(X,Y) are located at (I,J) elements
J=max(0,floor(n*(Y-eps)));
M=I+1+J*n; %(I,J) element is the M-th element

le=zeros(2,n^2); 
%   le: locating elements, two-row matrix. le(:,m) shows where the m-th
%   element locates. le(2,m)=1 (2, or 0) means that the m-th element is in
%   el1 (el2 or eli). le(1,m) is the corresponding column. That is,
%       if le(2,m)=1 then el1(le(1,m))=m,
%       if le(2,m)=2 then el2(le(1,m))=m,
%       if le(2,m)=0 then eli(le(1,m))=m,
le(2,el1)=1;
le(1,el1)=1:length(el1);
le(2,el2)=2;
le(1,el2)=1:length(el2);
le(1,eli(1,:))=1:length(eli(1,:));

%For each point on the interface, find an element that is not in Omega_2 and contains the point.
m=find(lvsfun(X,Y)==0);
k=find(le(2,M(m))==2);
I(m(k))=I(m(k))+1;
J(m(k))=J(m(k))+1;
M(m(k))=M(m(k))+n+1;

leM=le(:,M);
i1=find(leM(2,:)==1); %Elements in Omega_1
i2=find(leM(2,:)==2);
i0=find(leM(2,:)==0);
phi=lvsfun(X(i0),Y(i0));
i01=i0(phi<=0); %points in Omega_1 and on the interface
i02=i0(phi>0); %points in Omega_2

X=n*X-I; %Map (X,Y) into the reference element, still denoted by (X,Y).
Y=n*Y-J;

U=zeros(row*col,1);
for l=1:ldofs,
   jl=floor((l-1)/(p+1))+1;
   il=l-(jl-1)*(p+1);
   phi=upfemrefbas(X,Y,il,jl,p,xn);
   if ~isempty(i1),
       U(i1)=U(i1)+u(ltgmap1(l,leM(1,i1))).*phi(i1);
   end
   if ~isempty(i2),
       U(i2)=U(i2)+u(ltgmap2(l,leM(1,i2))).*phi(i2);
   end
   if ~isempty(i01),
       U(i01)=U(i01)+u(ltgmapi1(l,leM(1,i01))).*phi(i01);
   end
   if ~isempty(i02),
       U(i02)=U(i02)+u(ltgmapi2(l,leM(1,i02))).*phi(i02);   
   end
end
U=reshape(U,row,col);