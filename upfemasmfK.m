function [F1,F2,Fi1,Fi2]=upfemasmfK(a,f1,f2,gd,gn,lvsfun,lvsfunx,lvsfuny,n,p,beta,gamma0)
%Assemble element stiffness matrices
%   a: the coefficient in -div(a(x)grad u). a(x) is piecewise constant.
%       a=[a1 a2] where aj=a(x)|_{Omega_j}, j=1,2.
%   f1,f2: right hand sides in Omega_1 and Omega_2.
%   gd, gn: Drichlet and Neumann interface data.
%   lvsfun: its zero level set is the interface.
%   lvsfunx: the partial derivative of the lvsfun w.r.t. x.
%   lvsfuny: the partial derivative of the lvsfun w.r.t. y.
%   n: h=1/n is the mesh size
%   p: the degree of element polynomials (in each x and y, respectively).
%   Fj,j=1,2: local RHS vectors for elements in elj. 
%       The local DOFs in each element are sorted first along the x-axis then the y-axis.              
%   Fij,j=1,2: local RHS vectors (F_{hK}(phi_l^j))_{lm} for elements crossed 
%       by the interface. Fij(:,m) is the local RHS vector corresponding to the eli(m)-th element. 
%   xn: nodal points in [0,1].      
%       

h=1/n; %Mesh size
xn=linspace(0,1,p+1).';
ldofs=(p+1)^2; %Number of DOFs in reference element
w1=a(2)/(a(1)+a(2)); %Weights for average
w2=a(1)/(a(1)+a(2));
wa=a(1)*a(2)/(a(1)+a(2)); %w_1a_1=w_2a_2

%Default parameters for UPFEM
if nargin==10
    beta=1; %Symmetric
    gamma0=100;
end

%Find the elements in Omega_1, Omega_2, and the elements crossed by the interface.
[el1,el2,eli]=upfemeleide(lvsfun,lvsfunx,lvsfuny,n);

% Elements in Omega_1 and Omega_2
for k=1:2
    ks=num2str(k);
    elk=eval(['el' ks]);
    fk=eval(['f' ks]);
    eval(['F' ks '=zeros(ldofs,length(elk));']);
    jk=floor((elk-1)/n);
    ik=elk-1-jk*n;
    a=ik.'*h; 
    b=a+h;
    c=jk.'*h;
    d=c+h;
    for l=1:ldofs
        jl=floor((l-1)/(p+1))+1;
        il=l-(jl-1)*(p+1);
        f=@(x,y) fk(x,y).*upfemrefbas(x/h-floor(x/h),y/h-floor(y/h),il,jl,p,xn);
        eval(['F' ks '(l,:)=quad2(f,a,b,c,d,p);']);
    end
end
% Interface elements
el=eli(1,:);
jj=floor((el-1)/n);
ii=el-1-jj*n;
a=ii.'*h; 
b=a+h;
c=jj.'*h;
d=c+h;
Ni=size(eli,2);
Fi1=zeros(ldofs,Ni);
Fi2=zeros(ldofs,Ni);
for l=1:ldofs
    jl=floor((l-1)/(p+1))+1;
    il=l-(jl-1)*(p+1);
    %(fj,phi_l^j)_K
    f=@(x,y) f1(x,y).*upfemrefbas(x/h-floor(x/h),y/h-floor(y/h),il,jl,p,xn);
    Fi1(l,:)=upfemitfquad2d(f,lvsfun,lvsfunx,lvsfuny,eli,n,p,-1);
    f=@(x,y) f2(x,y).*upfemrefbas(x/h-floor(x/h),y/h-floor(y/h),il,jl,p,xn);
    Fi2(l,:)=upfemitfquad2d(f,lvsfun,lvsfunx,lvsfuny,eli,n,p,1);
    %<gn,{phi_l^j}^w>_e
    g=@(x,y) gn(x,y).*upfemrefbas(x/h-floor(x/h),y/h-floor(y/h),il,jl,p,xn);
    bl=upfemitfquad(g,lvsfun,lvsfunx,lvsfuny,eli,n,p);
    Fi1(l,:)=Fi1(l,:)+w2*bl;
    Fi2(l,:)=Fi2(l,:)+w1*bl;
    %beta<gd,{a*grad phi_j n}_w>_e
    g=@(x,y) gd(x,y)...
             .*(upfemrefbasx(x/h-floor(x/h),y/h-floor(y/h),il,jl,p,xn).*upfemitfnx(lvsfunx,lvsfuny,x,y)...
               +upfemrefbasy(x/h-floor(x/h),y/h-floor(y/h),il,jl,p,xn).*upfemitfny(lvsfunx,lvsfuny,x,y))/h;
    bbl=beta*upfemitfquad(g,lvsfun,lvsfunx,lvsfuny,eli,n,p);
    Fi1(l,:)=Fi1(l,:)-wa*bbl;
    Fi2(l,:)=Fi2(l,:)-wa*bbl;
    %gamma0{a}_w/h<gd,[phi_l^j]>_e
    g=@(x,y) gd(x,y).*upfemrefbas(x/h-floor(x/h),y/h-floor(y/h),il,jl,p,xn);
    j0l1=gamma0/h*2*wa*upfemitfquad(g,lvsfun,lvsfunx,lvsfuny,eli,n,p);
    Fi1(l,:)=Fi1(l,:)+j0l1;
    Fi2(l,:)=Fi2(l,:)-j0l1;
end
