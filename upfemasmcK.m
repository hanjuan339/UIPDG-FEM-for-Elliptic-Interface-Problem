function [K1,K2,Ki11,Ki12,Ki21,Ki22,MM,M12,M11,J12,J21,J11,J22]=upfemasmcK(a,theta,lvsfun,lvsfunx,lvsfuny,ell1,ell2,eli,e,n,p,beta,gamma0,gamma1)

h=1/n;
xn=linspace(0,1,p+1).';
ldofs=(p+1)^2;
wa=a(1)*a(2)/(a(1)+a(2));

if nargin==10
    beta=1;
    gamma0=100;
end

KK=zeros(ldofs,ldofs);
mm=zeros(ldofs,ldofs);
for k=1:ldofs
    jk=floor((k-1)/(p+1))+1;
    ik=k-(jk-1)*(p+1);
    for l=k:ldofs
        jl=floor((l-1)/(p+1))+1;
        il=l-(jl-1)*(p+1);
        KK(k,l)=quad2(@(x,y) upfemrefbasx(x,y,ik,jk,p,xn).*upfemrefbasx(x,y,il,jl,p,xn)...
                             +upfemrefbasy(x,y,ik,jk,p,xn).*upfemrefbasy(x,y,il,jl,p,xn),...
                              0,1,0,1,p);
      mm(k,l)=quad2(@(x,y) upfemrefbas(x,y,ik,jk,p,xn).*upfemrefbas(x,y,il,jl,p,xn)*h^2,0,1,0,1,p);
    end
end
KK=KK+triu(KK,1).';
mm=mm+triu(mm,1).';

K1=repmat(a(1)*KK(:)+theta(1)*mm(:),1,length(ell1));
K2=repmat(a(2)*KK(:)+theta(2)*mm(:),1,length(ell2));

A11=zeros((p+1)^4,length(eli)); %(grad phi_k^1, grad phi_l^1)_K1
                                %A22=(grad phi_k, grad phi_l)_K-A11. A12=A21=0.
B11=zeros((p+1)^4,length(eli)); %<{grad phi_k^1 cdot n}, [phi_l^1]>_e, 
                                %B12=-B11, B21=B11, B22=-B12
B11t=zeros((p+1)^4,length(eli)); %<[phi_k^1], {grad phi_l^1 cdot n}>_e,
                                %B12t=B11t,B21t=-B11t,B22t=-B11t
J011=zeros((p+1)^4,length(eli));%gamma0 {a}_w /h <[phi_k^1], [phi_l^1]>_e, 
                                %J012=-J011,J021=-J011,J022=J011
M11=zeros((p+1)^4,length(eli));

for k=1:ldofs
    jk=floor((k-1)/(p+1))+1;
    ik=k-(jk-1)*(p+1);
    for l=k:ldofs
        jl=floor((l-1)/(p+1))+1;
        il=l-(jl-1)*(p+1);
        %A11
        g=@(x,y) (upfemrefbasx(x/h-floor(x/h),y/h-floor(y/h),ik,jk,p,xn).*upfemrefbasx(x/h-floor(x/h),y/h-floor(y/h),il,jl,p,xn)...
                +upfemrefbasy(x/h-floor(x/h),y/h-floor(y/h),ik,jk,p,xn).*upfemrefbasy(x/h-floor(x/h),y/h-floor(y/h),il,jl,p,xn))/h^2;
        A11(k+ldofs*(l-1),:)=upfemitfquad2d(g,lvsfun,lvsfunx,lvsfuny,eli,n,p,-1);
        %M11
        g=@(x,y) upfemrefbas(x/h-floor(x/h),y/h-floor(y/h),ik,jk,p,xn).*upfemrefbas(x/h-floor(x/h),y/h-floor(y/h),il,jl,p,xn);
        M11(k+ldofs*(l-1),:)=upfemitfquad2d(g,lvsfun,lvsfunx,lvsfuny,eli,n,p,-1);
        %B11 and B11t
        g=@(x,y) (upfemrefbasx(x/h-floor(x/h),y/h-floor(y/h),ik,jk,p,xn).*upfemitfnx(lvsfunx,lvsfuny,x,y)...
                 +upfemrefbasy(x/h-floor(x/h),y/h-floor(y/h),ik,jk,p,xn).*upfemitfny(lvsfunx,lvsfuny,x,y))...
                 .*upfemrefbas(x/h-floor(x/h),y/h-floor(y/h),il,jl,p,xn)/h;
        B11(k+ldofs*(l-1),:)=upfemitfquad(g,lvsfun,lvsfunx,lvsfuny,eli,n,p);
        B11t(l+ldofs*(k-1),:)=B11(k+ldofs*(l-1),:);
        %J011
        g=@(x,y) upfemrefbas(x/h-floor(x/h),y/h-floor(y/h),ik,jk,p,xn).*upfemrefbas(x/h-floor(x/h),y/h-floor(y/h),il,jl,p,xn);
        J011(k+ldofs*(l-1),:)=gamma0/h*2*wa*upfemitfquad(g,lvsfun,lvsfunx,lvsfuny,eli,n,p);
        if l~=k
            A11(l+ldofs*(k-1),:)=A11(k+ldofs*(l-1),:);
            g=@(x,y) (upfemrefbasx(x/h-floor(x/h),y/h-floor(y/h),il,jl,p,xn).*upfemitfnx(lvsfunx,lvsfuny,x,y)...
                     +upfemrefbasy(x/h-floor(x/h),y/h-floor(y/h),il,jl,p,xn).*upfemitfny(lvsfunx,lvsfuny,x,y))...
                     .*upfemrefbas(x/h-floor(x/h),y/h-floor(y/h),ik,jk,p,xn)/h;
            B11(l+ldofs*(k-1),:)=upfemitfquad(g,lvsfun,lvsfunx,lvsfuny,eli,n,p);
            B11t(k+ldofs*(l-1),:)=B11(l+ldofs*(k-1),:);
            J011(l+ldofs*(k-1),:)=J011(k+ldofs*(l-1),:);
        end
    end
end
Ki11=a(1)*A11+theta(1)*M11+J011-wa*(B11+beta*B11t);
Ki22=a(2)*(repmat(KK(:),1,length(eli))-A11)+theta(2)*(repmat(mm(:),1,length(eli))-M11)+J011+wa*(B11+beta*B11t);
Ki12=-J011+wa*(B11-beta*B11t);
Ki21=-J011-wa*(B11-beta*B11t);



% 整边加罚项
MM=zeros(p+1);
for k=1:p+1
    for l=1:p+1
        MM(k,l)=quad1(@(x,s) upfemrefbas(x,0*x,k,1,p,xn).*upfemrefbas(x,0*x,l,1,p,xn),0,1,p);
    end
end
MM=MM*gamma1;%这里det(B)产生的h和格式里的1/h抵消了
    
% 整边 适定项
M12=cell(1,4);
M12{1}=zeros(ldofs);
M12{2}=M12{1};
M12{3}=M12{1};
M12{4}=M12{1};



for k=1:ldofs
    jk=floor((k-1)/(p+1))+1;
    ik=k-(jk-1)*(p+1);
    for l=1:ldofs
        jl=floor((l-1)/(p+1))+1;
        il=l-(jl-1)*(p+1);
        M12{1}(k,l)=quad1(@(x,s) upfemrefbasy(x,0*x,ik,jk,p,xn).*upfemrefbas(x,0*x+1,il,jl,p,xn)/2 ...
            -beta*upfemrefbas(x,0*x,ik,jk,p,xn).*upfemrefbasy(x,0*x+1,il,jl,p,xn)/2,0,1,p);
        M12{2}(k,l)=quad1(@(y,s) -upfemrefbasx(0*y+1,y,ik,jk,p,xn).*upfemrefbas(0*y,y,il,jl,p,xn)/2 ...
            +beta*upfemrefbas(0*y+1,y,ik,jk,p,xn).*upfemrefbasx(0*y,y,il,jl,p,xn)/2,0,1,p);
        M12{3}(k,l)=quad1(@(x,s) -upfemrefbasy(x,0*x+1,ik,jk,p,xn).*upfemrefbas(x,0*x,il,jl,p,xn)/2 ...
            +beta*upfemrefbas(x,0*x+1,ik,jk,p,xn).*upfemrefbasy(x,0*x,il,jl,p,xn)/2,0,1,p);
        M12{4}(k,l)=quad1(@(y,s) upfemrefbasx(0*y,y,ik,jk,p,xn).*upfemrefbas(0*y+1,y,il,jl,p,xn)/2 ...
            -beta*upfemrefbas(0*y,y,ik,jk,p,xn).*upfemrefbasx(0*y+1,y,il,jl,p,xn)/2,0,1,p);
    end
end

% for k=1:ldofs
%     jk=floor((k-1)/(p+1))+1;
%     ik=k-(jk-1)*(p+1);
%     for l=1:ldofs
%         jl=floor((l-1)/(p+1))+1;
%         il=l-(jl-1)*(p+1);
%         M21{1}(k,l)=quad1(@(x) -upfemrefbasy(x,0*x+1,ik,jk,p,xn).*upfemrefbas(x,0*x,il,jl,p,xn)/2 ...
%             +beta*upfemrefbas(x,0*x+1,ik,jk,p,xn).*upfemrefbasy(x,0*x,il,jl,p,xn)/2,0,1,p);
%         M21{2}(k,l)=quad1(@(y) upfemrefbasx(0*y,y,ik,jk,p,xn).*upfemrefbas(0*y+1,y,il,jl,p,xn)/2 ...
%             -beta*upfemrefbas(0*y,y,ik,jk,p,xn).*upfemrefbasx(0*y+1,y,il,jl,p,xn)/2,0,1,p);
%         M21{3}(k,l)=quad1(@(x) upfemrefbasy(x,0*x,ik,jk,p,xn).*upfemrefbas(x,0*x+1,il,jl,p,xn)/2 ...
%             -beta*upfemrefbas(x,0*x,ik,jk,p,xn).*upfemrefbasy(x,0*x+1,il,jl,p,xn)/2,0,1,p);
%         M21{4}(k,l)=quad1(@(y) -upfemrefbasx(0*y+1,y,ik,jk,p,xn).*upfemrefbas(0*y,y,il,jl,p,xn)/2 ...
%             +beta*upfemrefbas(0*y+1,y,ik,jk,p,xn).*upfemrefbasx(0*y,y,il,jl,p,xn)/2,0,1,p);
%     end
% end

M11=cell(1,4);
M11{1}=zeros(ldofs);
M11{2}=M11{1};
M11{3}=M11{1};
M11{4}=M11{1};

for k=1:ldofs
    jk=floor((k-1)/(p+1))+1;
    ik=k-(jk-1)*(p+1);
    for l=1:ldofs
        jl=floor((l-1)/(p+1))+1;
        il=l-(jl-1)*(p+1);
        M11{1}(k,l)=quad1(@(x,s) -upfemrefbasy(x,0*x,ik,jk,p,xn).*upfemrefbas(x,0*x,il,jl,p,xn)/2 ...
            -beta*upfemrefbas(x,0*x,ik,jk,p,xn).*upfemrefbasy(x,0*x,il,jl,p,xn)/2,0,1,p);
        M11{2}(k,l)=quad1(@(y,s) upfemrefbasx(0*y+1,y,ik,jk,p,xn).*upfemrefbas(0*y+1,y,il,jl,p,xn)/2 ...
            +beta*upfemrefbas(0*y+1,y,ik,jk,p,xn).*upfemrefbasx(0*y+1,y,il,jl,p,xn)/2,0,1,p);
        M11{3}(k,l)=quad1(@(x,s) upfemrefbasy(x,0*x+1,ik,jk,p,xn).*upfemrefbas(x,0*x+1,il,jl,p,xn)/2 ...
            +beta*upfemrefbas(x,0*x+1,ik,jk,p,xn).*upfemrefbasy(x,0*x+1,il,jl,p,xn)/2,0,1,p);
        M11{4}(k,l)=quad1(@(y,s) -upfemrefbasx(0*y,y,ik,jk,p,xn).*upfemrefbas(0*y,y,il,jl,p,xn)/2 ...
            -beta*upfemrefbas(0*y,y,ik,jk,p,xn).*upfemrefbasx(0*y,y,il,jl,p,xn)/2,0,1,p);
    end
end



% 断边加罚项和适定项


index=[-1 1];

J12=cell(1,2);
J21=cell(1,2);
J11=cell(1,2);
J22=cell(1,2);
for l=1:2
    ee=e{2*l};
    y1=floor((ee(1,:)-1)/n);
    x1=ee(1,:)-n*y1-1;
    x2=x1+1;y2=y1;
    x3=x2;y3=y2+1;
    x4=x1;y4=y1+1;
    d=[lvsfun(x1*h,y1*h);lvsfun(x2*h,y2*h);lvsfun(x3*h,y3*h);lvsfun(x4*h,y4*h)]*index(l);
    J12{l}=zeros((p+1)^4,4*size(ee,2));
    J21{l}=zeros((p+1)^4,4*size(ee,2));
    J11{l}=zeros((p+1)^4,4*size(ee,2));
    J22{l}=zeros((p+1)^4,4*size(ee,2));
    x1=x1*h;y1=y1*h;
    for j=1:size(ee,2)
        reflvsfun=@(x,y) lvsfun(x1(j)+h*x,y1(j)+h*y);
        reflvsfunx=@(x,y) h*lvsfunx(x1(j)+h*x,y1(j)+h*y);
        reflvsfuny=@(x,y) h*lvsfuny(x1(j)+h*x,y1(j)+h*y);
        %i=1
        if ee(2,j)~=0
            for k=1:(p+1)^2
                kj=floor((k-1)/(p+1))+1;
                ki=k-(p+1)*(kj-1);
                for m=1:(p+1)^2
                    mj=floor((m-1)/(p+1))+1;
                    mi=m-(p+1)*(mj-1);
                    %J12
                    g12=@(x,s) -upfemrefbas(x,0*x,ki,kj,p,xn).*upfemrefbas(x,0*x+1,mi,mj,p,xn)*gamma1*a(l)...
                        -upfemrefbasy(x,0*x,ki,kj,p,xn).*upfemrefbas(x,0*x+1,mi,mj,p,xn)*a(l)/2 ...
                        +beta*upfemrefbas(x,0*x,ki,kj,p,xn).*upfemrefbasy(x,0*x+1,mi,mj,p,xn)*a(l)/2;
                    %J21
                    g21=@(x,s) -upfemrefbas(x,0*x+1,ki,kj,p,xn).*upfemrefbas(x,0*x,mi,mj,p,xn)*gamma1*a(l)...
                        +upfemrefbasy(x,0*x+1,ki,kj,p,xn).*upfemrefbas(x,0*x,mi,mj,p,xn)*a(l)/2 ...
                        -beta*upfemrefbas(x,0*x+1,ki,kj,p,xn).*upfemrefbasy(x,0*x,mi,mj,p,xn)*a(l)/2;
                    %J11
                    g11=@(x,s) upfemrefbas(x,0*x,ki,kj,p,xn).*upfemrefbas(x,0*x,mi,mj,p,xn)*gamma1*a(l)...
                        +upfemrefbasy(x,0*x,ki,kj,p,xn).*upfemrefbas(x,0*x,mi,mj,p,xn)*a(l)/2 ...
                        +beta*upfemrefbas(x,0*x,ki,kj,p,xn).*upfemrefbasy(x,0*x,mi,mj,p,xn)*a(l)/2;
                    %J22
                    g22=@(x,s) upfemrefbas(x,0*x+1,ki,kj,p,xn).*upfemrefbas(x,0*x+1,mi,mj,p,xn)*gamma1*a(l)...
                        -upfemrefbasy(x,0*x+1,ki,kj,p,xn).*upfemrefbas(x,0*x+1,mi,mj,p,xn)*a(l)/2 ...
                        -beta*upfemrefbas(x,0*x+1,ki,kj,p,xn).*upfemrefbasy(x,0*x+1,mi,mj,p,xn)*a(l)/2;
                    
                    x0=fzerox(reflvsfun,reflvsfunx,0.5,0);
                    if d(1,j)>0
                        J12{l}((m-1)*(p+1)^2+k,4*(j-1)+1)=quad1(g12,0,x0,p);
                        J21{l}((m-1)*(p+1)^2+k,4*(j-1)+1)=quad1(g21,0,x0,p);
                        J11{l}((m-1)*(p+1)^2+k,4*(j-1)+1)=quad1(g11,0,x0,p);
                        J22{l}((m-1)*(p+1)^2+k,4*(j-1)+1)=quad1(g22,0,x0,p);
                    else
                        J12{l}((m-1)*(p+1)^2+k,4*(j-1)+1)=quad1(g12,x0,1,p);
                        J21{l}((m-1)*(p+1)^2+k,4*(j-1)+1)=quad1(g21,x0,1,p);
                        J11{l}((m-1)*(p+1)^2+k,4*(j-1)+1)=quad1(g11,x0,1,p);
                        J22{l}((m-1)*(p+1)^2+k,4*(j-1)+1)=quad1(g22,x0,1,p);
                    end
                   
                end
            end
        end
        %i=2
        if ee(3,j)~=0
            for k=1:(p+1)^2
                kj=floor((k-1)/(p+1))+1;
                ki=k-(p+1)*(kj-1);
                for m=1:(p+1)^2
                    mj=floor((m-1)/(p+1))+1;
                    mi=m-(p+1)*(mj-1);
                    %J12
                    g12=@(y,s) -upfemrefbas(0*y+1,y,ki,kj,p,xn).*upfemrefbas(0*y,y,mi,mj,p,xn)*gamma1*a(l)...
                        +upfemrefbasx(0*y+1,y,ki,kj,p,xn).*upfemrefbas(0*y,y,mi,mj,p,xn)*a(l)/2 ...
                        -beta*upfemrefbas(0*y+1,y,ki,kj,p,xn).*upfemrefbasx(0*y,y,mi,mj,p,xn)*a(l)/2;
                    g21=@(y,s) -upfemrefbas(0*y,y,ki,kj,p,xn).*upfemrefbas(0*y+1,y,mi,mj,p,xn)*gamma1*a(l)...
                        -upfemrefbasx(0*y,y,ki,kj,p,xn).*upfemrefbas(0*y+1,y,mi,mj,p,xn)*a(l)/2 ...
                        +beta*upfemrefbas(0*y,y,ki,kj,p,xn).*upfemrefbasx(0*y+1,y,mi,mj,p,xn)*a(l)/2;
                    g11=@(y,s) upfemrefbas(0*y+1,y,ki,kj,p,xn).*upfemrefbas(0*y+1,y,mi,mj,p,xn)*gamma1*a(l)...
                        -upfemrefbasx(0*y+1,y,ki,kj,p,xn).*upfemrefbas(0*y+1,y,mi,mj,p,xn)*a(l)/2 ...
                        -beta*upfemrefbas(0*y+1,y,ki,kj,p,xn).*upfemrefbasx(0*y+1,y,mi,mj,p,xn)*a(l)/2;
                    g22=@(y,s) upfemrefbas(0*y,y,ki,kj,p,xn).*upfemrefbas(0*y,y,mi,mj,p,xn)*gamma1*a(l)...
                        +upfemrefbasx(0*y,y,ki,kj,p,xn).*upfemrefbas(0*y,y,mi,mj,p,xn)*a(l)/2 ...
                        +beta*upfemrefbas(0*y,y,ki,kj,p,xn).*upfemrefbasx(0*y,y,mi,mj,p,xn)*a(l)/2;
                    x0=fzeroy(reflvsfun,reflvsfuny,1,0.5);
                    if d(2,j)>0
                        J12{l}((m-1)*(p+1)^2+k,4*(j-1)+2)=quad1(g12,0,x0,p);
                        J21{l}((m-1)*(p+1)^2+k,4*(j-1)+2)=quad1(g21,0,x0,p);
                        J11{l}((m-1)*(p+1)^2+k,4*(j-1)+2)=quad1(g11,0,x0,p);
                        J22{l}((m-1)*(p+1)^2+k,4*(j-1)+2)=quad1(g22,0,x0,p);
                    else
                        J12{l}((m-1)*(p+1)^2+k,4*(j-1)+2)=quad1(g12,x0,1,p);
                        J21{l}((m-1)*(p+1)^2+k,4*(j-1)+2)=quad1(g21,x0,1,p);
                        J11{l}((m-1)*(p+1)^2+k,4*(j-1)+2)=quad1(g11,x0,1,p);
                        J22{l}((m-1)*(p+1)^2+k,4*(j-1)+2)=quad1(g22,x0,1,p);
                    end
                   
                end
            end
        end
        %i=3
        if ee(4,j)~=0
            for k=1:(p+1)^2
                kj=floor((k-1)/(p+1))+1;
                ki=k-(p+1)*(kj-1);
                for m=1:(p+1)^2
                    mj=floor((m-1)/(p+1))+1;
                    mi=m-(p+1)*(mj-1);
                    %J12
                    g12=@(x,s) -upfemrefbas(x,0*x+1,ki,kj,p,xn).*upfemrefbas(x,0*x,mi,mj,p,xn)*gamma1*a(l)...
                        +upfemrefbasy(x,0*x+1,ki,kj,p,xn).*upfemrefbas(x,0*x,mi,mj,p,xn)*a(l)/2 ...
                        -beta*upfemrefbas(x,0*x+1,ki,kj,p,xn).*upfemrefbasy(x,0*x,mi,mj,p,xn)*a(l)/2;
                    g21=@(x,s) -upfemrefbas(x,0*x,ki,kj,p,xn).*upfemrefbas(x,0*x+1,mi,mj,p,xn)*gamma1*a(l)...
                        -upfemrefbasy(x,0*x,ki,kj,p,xn).*upfemrefbas(x,0*x+1,mi,mj,p,xn)*a(l)/2 ...
                        +beta*upfemrefbas(x,0*x,ki,kj,p,xn).*upfemrefbasy(x,0*x+1,mi,mj,p,xn)*a(l)/2;
                    g11=@(x,s) upfemrefbas(x,0*x+1,ki,kj,p,xn).*upfemrefbas(x,0*x+1,mi,mj,p,xn)*gamma1*a(l)...
                        -upfemrefbasy(x,0*x+1,ki,kj,p,xn).*upfemrefbas(x,0*x+1,mi,mj,p,xn)*a(l)/2 ...
                        -beta*upfemrefbas(x,0*x+1,ki,kj,p,xn).*upfemrefbasy(x,0*x+1,mi,mj,p,xn)*a(l)/2;
                    g22=@(x,s) upfemrefbas(x,0*x,ki,kj,p,xn).*upfemrefbas(x,0*x,mi,mj,p,xn)*gamma1*a(l)...
                        +upfemrefbasy(x,0*x,ki,kj,p,xn).*upfemrefbas(x,0*x,mi,mj,p,xn)*a(l)/2 ...
                        +beta*upfemrefbas(x,0*x,ki,kj,p,xn).*upfemrefbasy(x,0*x,mi,mj,p,xn)*a(l)/2;
                    x0=fzerox(reflvsfun,reflvsfunx,0.5,1);
                    if d(4,j)>0
                        J12{l}((m-1)*(p+1)^2+k,4*(j-1)+3)=quad1(g12,0,x0,p);
                        J21{l}((m-1)*(p+1)^2+k,4*(j-1)+3)=quad1(g21,0,x0,p);
                        J11{l}((m-1)*(p+1)^2+k,4*(j-1)+3)=quad1(g11,0,x0,p);
                        J22{l}((m-1)*(p+1)^2+k,4*(j-1)+3)=quad1(g22,0,x0,p);
                    else
                        J12{l}((m-1)*(p+1)^2+k,4*(j-1)+3)=quad1(g12,x0,1,p);
                        J21{l}((m-1)*(p+1)^2+k,4*(j-1)+3)=quad1(g21,x0,1,p);
                        J11{l}((m-1)*(p+1)^2+k,4*(j-1)+3)=quad1(g11,x0,1,p);
                        J22{l}((m-1)*(p+1)^2+k,4*(j-1)+3)=quad1(g22,x0,1,p);
                    end
                    
                end
            end
        end
        %i=4
        if ee(5,j)~=0
            for k=1:(p+1)^2
                kj=floor((k-1)/(p+1))+1;
                ki=k-(p+1)*(kj-1);
                for m=1:(p+1)^2
                    mj=floor((m-1)/(p+1))+1;
                    mi=m-(p+1)*(mj-1);
                    %J12
                    g12=@(y,s) -upfemrefbas(0*y,y,ki,kj,p,xn).*upfemrefbas(0*y+1,y,mi,mj,p,xn)*gamma1*a(l)...
                        -upfemrefbasx(0*y,y,ki,kj,p,xn).*upfemrefbas(0*y+1,y,mi,mj,p,xn)*a(l)/2 ...
                        +beta*upfemrefbas(0*y,y,ki,kj,p,xn).*upfemrefbasx(0*y+1,y,mi,mj,p,xn)*a(l)/2;
                    g21=@(y,s) -upfemrefbas(0*y+1,y,ki,kj,p,xn).*upfemrefbas(0*y,y,mi,mj,p,xn)*gamma1*a(l)...
                        +upfemrefbasx(0*y+1,y,ki,kj,p,xn).*upfemrefbas(0*y,y,mi,mj,p,xn)*a(l)/2 ...
                        -beta*upfemrefbas(0*y+1,y,ki,kj,p,xn).*upfemrefbasx(0*y,y,mi,mj,p,xn)*a(l)/2;
                    g11=@(y,s) upfemrefbas(0*y,y,ki,kj,p,xn).*upfemrefbas(0*y,y,mi,mj,p,xn)*gamma1*a(l)...
                        +upfemrefbasx(0*y,y,ki,kj,p,xn).*upfemrefbas(0*y,y,mi,mj,p,xn)*a(l)/2 ...
                        +beta*upfemrefbas(0*y,y,ki,kj,p,xn).*upfemrefbasx(0*y,y,mi,mj,p,xn)*a(l)/2;
                    g22=@(y,s) upfemrefbas(0*y+1,y,ki,kj,p,xn).*upfemrefbas(0*y+1,y,mi,mj,p,xn)*gamma1*a(l)...
                        -upfemrefbasx(0*y+1,y,ki,kj,p,xn).*upfemrefbas(0*y+1,y,mi,mj,p,xn)*a(l)/2 ...
                        -beta*upfemrefbas(0*y+1,y,ki,kj,p,xn).*upfemrefbasx(0*y+1,y,mi,mj,p,xn)*a(l)/2;
                    x0=fzeroy(reflvsfun,reflvsfuny,0,0.5);
                    if d(1,j)>0
                        J12{l}((m-1)*(p+1)^2+k,4*(j-1)+4)=quad1(g12,0,x0,p);
                        J21{l}((m-1)*(p+1)^2+k,4*(j-1)+4)=quad1(g21,0,x0,p);
                        J11{l}((m-1)*(p+1)^2+k,4*(j-1)+4)=quad1(g11,0,x0,p);
                        J22{l}((m-1)*(p+1)^2+k,4*(j-1)+4)=quad1(g22,0,x0,p);
                    else
                        J12{l}((m-1)*(p+1)^2+k,4*(j-1)+4)=quad1(g12,x0,1,p);
                        J21{l}((m-1)*(p+1)^2+k,4*(j-1)+4)=quad1(g21,x0,1,p);
                        J11{l}((m-1)*(p+1)^2+k,4*(j-1)+4)=quad1(g11,x0,1,p);
                        J22{l}((m-1)*(p+1)^2+k,4*(j-1)+4)=quad1(g22,x0,1,p);
                    end
                    
                end
            end
        end
    end
end

