function U=upfemeleevaluh(x,y,u,ltgmap,n,p,xn)
%Evaluating the UPFE soltion on given points
%   x,y: the coordinates of the points. x,y may be matrics of the same size.
%   ltgmap:  ltgmap(j,m) is the index to the golbal DOF of the j-th local
%       DOF in the m-th element.         
%   u: Nodal values of the discrete solution.
%   xn: nodal points in [0,1].      

ldofs=(p+1)^2;
[row,col]=size(x);
X=x(:);
Y=y(:);
I=max(0,floor(n*(X-eps))); %(X,Y) are located at (I,J) elements
J=max(0,floor(n*(Y-eps)));
el=I+1+n*J; % Element index
ltgmap=ltgmap(:,el);
X=n*X-I; %Map (X,Y) into the reference element, still denoted by (X,Y).
Y=n*Y-J;

U=zeros(row*col,1);
for l=1:ldofs
   jl=floor((l-1)/(p+1))+1;
   il=l-(jl-1)*(p+1);
   phi=upfemrefbas(X,Y,il,jl,p,xn);
   U=U+u(ltgmap(l,:)).*phi;
end
U=reshape(U,row,col);