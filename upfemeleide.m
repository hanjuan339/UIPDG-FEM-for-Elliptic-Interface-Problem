function [el1,el2,eli]=upfemeleide(lvsfun,lvsfunx,lvsfuny,n)
%Identify elements
%   lvsfun: specifies the interface (zero level set). It must be
%   vectorized.
%   n: The unit square domain is divided into n by n sub-squares
%   el1: elements inside the interface (where the lvsfun<0)
%   el2: elements outside the interface (where the lvsfun>0)
%   eli: interface elements. eli is a six-row matrix. The first row gives
%       the indices. The second row indicates typies of elements: 
%        1: The interface cuts the left and right edges, the upper edge in Omega_1;
%        2: cuts the left and right edges, the lower edge in Omega_1;
%        3: cuts the lower and upper edges, the right edge in Omega_1;
%        4: cuts the lower and upper edges, the left edge in Omega_1;
%        5: cuts the upper and right edges, the upper-right conner in Omega_1;
%        6: cuts the upper and left edges, the upper-left conner in Omega_1;
%        7: cuts the lower and right edges, the lower-right conner in Omega_1;                  
%        8: cuts the lower and left edges, the lower-left conner in Omega_1;
%        9: cuts the lower and left edges, the lower-left conner in Omega_2;
%        10: cuts the lower and right edges, the lower-right conner in Omega_2;
%        11: cuts the upper and left edges, the upper-left conner in Omega_2;
%        12: cuts the upper and right edges, the upper-right conner in Omega_2;
%        13: cuts the upper edge twice, four conners in Omega_2;
%        14: cuts the lower edge twice and the left and right edges once, respectively, the upper edge in Omega_1;
%        15: cuts the lower edge twice, four conners in Omega_2;
%        16: cuts the upper edge twice and the left and right edges once, respectively, the lower edge in Omega_1;
%        17: cuts the right edge twice, four conners in Omega_2;
%        18: cuts the left edge twice and the lower and upper edges once, respectively, the right edge in Omega_1;
%        19: cuts the left edge twice, four conners in Omega_2;
%        20: cuts the right edge twice and the lower and upper edges once, respectively, the left edge in Omega_1;
%        21: cuts the upper edge twice, four conners in Omega_1;
%        22: cuts the lower edge twice and the left and right edges once, respectively, the upper edge in Omega_2;
%        23: cuts the lower edge twice, four conners in Omega_1;
%        24: cuts the upper edge twice and the left and right edges once, respectively, the lower edge in Omega_2;
%        25: cuts the right edge twice, four conners in Omega_1;
%        26: cuts the left edge twice and the lower and upper edges once, respectively, the right edge in Omega_2;
%        27: cuts the left edge twice, four conners in Omega_1;
%        28: cuts the right edge twice and the lower and upper edges once, respectively, the left edge in Omega_2.
%       The 3rd-4th rows stores the x-coordinates of the points of intersection (if exist) at the lower and upper edges;
%        or the x-coordinates of left and right points of intersection for type 13-18, or 21-24 elements.
%       The 5-6th rows stores the y-coordinates of the points of intersection (if exist) at the left and right edges;
%        or the y-coordinates of lower and upper points of intersection for type 17-20, or 25-28 elements.
%   The elements are ordered first along the x-axis then the y-axis 

h=1/n;

elitype=zeros(1,n^2);
eliedge=zeros(4,n^2); % The j-th column indicates the points of intersection 
                      % of the interface and the four edges of the j-th element.
                      % The edges are ordered as [lower;upper;left;right].
                      % For type 13-18, or 21-24 element: The left and right points of intersection on
                      %     the upper or lower edge are stored in the 1st and 2nd rows respectively.
                      % For type 17-20, or 25-28 element: The lower and upper points of intersection on
                      %     the left or right edge are stored in the 3rd and 4th rows respectively.
c=linspace(0,1,n+1);
[x,y]=meshgrid(c);
z=lvsfun(x,y);
[J,I]=find(z(1:n,:).*z(2:n+1,:)<=0);
% J=n+1-J;
y1=fzeroy(lvsfun,lvsfuny,(I-1)*h,J*h-h/2);
eliedge(4,(J-1)*n+I-1)=y1;
eliedge(3,(J-1)*n+I)=y1;
A=sortrows([I J y1],2);
i=find(A(2:end,2)-A(1:end-1,2)==0 & A(2:end,1)-A(1:end-1,1)==1);
nx1=lvsfunx((A(i,1)-1)*h,A(i,3));
nx2=lvsfunx((A(i+1,1)-1)*h,A(i+1,3));
ny1=lvsfuny((A(i,1)-1)*h,A(i,3));
ny2=lvsfuny((A(i+1,1)-1)*h,A(i+1,3));
j=find(ny1.*ny2>eps & nx1.*nx2<-eps);
for k=1:length(j)
    ijk=i(j(k));
    Ik=A(ijk,1);
    Jk=A(ijk,2);
    f=@(x) lvsfunx(x,fzeroy(lvsfun,lvsfuny,x,A(ijk,3)));
    xc=fzero(f,[(Ik-1)*h Ik*h]);
    yc=fzeroy(lvsfun,lvsfuny,xc,A(ijk,3));
    if  yc<(Jk-1)*h
        el=(Jk-2)*n+Ik;
        if ny1(j(k))<-eps
            elitype(el)=13;
            elitype(el+n)=14;
        else
            elitype(el)=21;
            elitype(el+n)=22;
        end
        xl=fzerox(lvsfun,lvsfunx,(Ik-1)*h,(Jk-1)*h);
        xr=fzerox(lvsfun,lvsfunx,Ik*h,(Jk-1)*h);
        eliedge(1,el)=xl;
        eliedge(2,el)=xr;
        eliedge(1,el+n)=xl;
        eliedge(2,el+n)=xr;
    end
    if  yc>Jk*h
        el=(Jk-1)*n+Ik;
        if ny1(j(k))>eps
            elitype(el)=16;
            elitype(el+n)=15;
        else
            elitype(el)=24;
            elitype(el+n)=23;
        end
        xl=fzerox(lvsfun,lvsfunx,(Ik-1)*h,Jk*h);
        xr=fzerox(lvsfun,lvsfunx,Ik*h,Jk*h);
        eliedge(1,el)=xl;
        eliedge(2,el)=xr;
        eliedge(1,el+n)=xl;
        eliedge(2,el+n)=xr;
    end
end
[J,I]=find(z(:,1:n).*z(:,2:n+1)<=0);
% J=n+2-J;
x1=fzerox(lvsfun,lvsfunx,I*h-h/2,(J-1)*h);
eliedge(2,(J-2)*n+I)=x1;
eliedge(1,(J-1)*n+I)=x1;
A=sortrows([I J x1],1);
i=find(A(1:end-1,1)-A(2:end,1)==0 & A(1:end-1,2)-A(2:end,2)==-1);
nx1=lvsfunx(A(i,3),(A(i,2)-1)*h);
nx2=lvsfunx(A(i+1,3),(A(i+1,2)-1)*h);
ny1=lvsfuny(A(i,3),(A(i,2)-1)*h);
ny2=lvsfuny(A(i+1,3),(A(i+1,2)-1)*h);
j=find(nx1.*nx2>eps & ny1.*ny2<-eps);
for k=1:length(j)
    ijk=i(j(k));
    Ik=A(ijk,1);
    Jk=A(ijk,2);
    f=@(y) lvsfuny(fzerox(lvsfun,lvsfunx,A(ijk,3),y),y);
    yc=fzero(f,[(Jk-1)*h Jk*h]);
    xc=fzerox(lvsfun,lvsfunx,A(ijk,3),yc);
    if  xc<(Ik-1)*h
        el=(Jk-1)*n+Ik-1;
        if nx1(j(k))<-eps
            elitype(el)=17;
            elitype(el+1)=18;
        else
            elitype(el)=25;
            elitype(el+1)=26;
        end
        yl=fzeroy(lvsfun,lvsfuny,(Ik-1)*h,(Jk-1)*h);
        yu=fzeroy(lvsfun,lvsfuny,(Ik-1)*h,Jk*h);
        eliedge(3,el)=yl;
        eliedge(4,el)=yu;
        eliedge(3,el+1)=yl;
        eliedge(4,el+1)=yu;
    end
    if  xc>Ik*h
        el=(Jk-1)*n+Ik;
        if nx1(j(k))>eps
            elitype(el)=20;
            elitype(el+1)=19;
        else
            elitype(el)=28;
            elitype(el+1)=27;
        end
        yl=fzeroy(lvsfun,lvsfuny,Ik*h,(Jk-1)*h);
        yu=fzeroy(lvsfun,lvsfuny,Ik*h,Jk*h);
        eliedge(3,el)=yl;
        eliedge(4,el)=yu;
        eliedge(3,el+1)=yl;
        eliedge(4,el+1)=yu;
    end
end

% Elements cutted by the interface
m=find(sum(eliedge~=0)>0);
Jm=floor((m-1)/n)+1;
Im=m-(Jm-1)*n;
eli=[m;elitype(m);eliedge(:,m)];
i=find((eli(3,:)==0 | eli(4,:)==0) & eli(5,:)~=0 & eli(6,:)~=0 & eli(2,:)==0);
j=z((Im(i)-1)*(n+1)+Jm(i))>-eps & z(Im(i)*(n+1)+Jm(i))>-eps;
eli(2,i(j)) = 1;
eli(2,i(~j)) = 2;
i=find(eli(3,:)~=0 & eli(4,:)~=0 & (eli(5,:)==0 | eli(6,:)==0) & eli(2,:)==0);
j=z((Im(i)-1)*(n+1)+Jm(i))>-eps & z((Im(i)-1)*(n+1)+Jm(i)+1)>-eps;
eli(2,i(j)) = 3;
eli(2,i(~j)) = 4;
i=find(eli(3,:)==0 & eli(4,:)~=0 & eli(5,:)==0 & eli(6,:)~=0 & eli(2,:)==0);
j=z((Im(i))*(n+1)+Jm(i)+1)<-eps;
eli(2,i(j)) = 5;
eli(2,i(~j)) = 12;
i=find(eli(3,:)==0 & eli(4,:)~=0 & eli(5,:)~=0 & eli(6,:)==0 & eli(2,:)==0);
j=z((Im(i)-1)*(n+1)+Jm(i)+1)<-eps;
eli(2,i(j)) = 6;
eli(2,i(~j)) = 11;
i=find(eli(3,:)~=0 & eli(4,:)==0 & eli(5,:)==0 & eli(6,:)~=0 & eli(2,:)==0);
j=z(Im(i)*(n+1)+Jm(i))<-eps;
eli(2,i(j)) = 7;
eli(2,i(~j)) = 10;
i=find(eli(3,:)~=0 & eli(4,:)==0 & eli(5,:)~=0 & eli(6,:)==0 & eli(2,:)==0);
j=z((Im(i)-1)*(n+1)+Jm(i))<-eps;
eli(2,i(j)) = 8;
eli(2,i(~j)) = 9;

[J,I]=find(z(1:n,1:n)<eps & z(1:n,2:n+1)<eps & z(2:n+1,1:n)<eps & z(2:n+1,2:n+1)<eps);
% J=n+1-J;
el=zeros(1,n^2);
el((J-1)*n+I)=-1; % Mark elements in Omega_1
[J,I]=find(z(1:n,1:n)>-eps & z(1:n,2:n+1)>-eps & z(2:n+1,1:n)>-eps & z(2:n+1,2:n+1)>-eps);
% J=n+1-J;
el((J-1)*n+I)=1; % Mark elements in Omega_2
el(m(eli(2,:)>12))=0; % Remove those interface elements whose four vertices in the same subdomain
el1=find(el<0);
el2=find(el>0);
eli(:,el(m)~=0)=[]; % Remove elements with only one vertex on the interface
