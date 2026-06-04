function q=upfemitfquad(g,lvsfun,lvsfunx,lvsfuny,eli,n,p)
%Evaluate line integral along an interface segment cotained in the unit suqare.
%   g: the integrand, a function of x and y.
%   lvsfun: the level set function whose zero level set defines the interface.
%   lvsfunx: the partial derivative of the lvsfun w.r.t. x
%   lvsfuny: the partial derivative of the lvsfun w.r.t. y
%   eli: Interface elements.
%   q: the integral of g.

ne=size(eli,2); % Number of interface elements
h=1/n;

% Conners of interface elements
el=eli(1,:);
jl=floor((el-1)/n);
il=el-1-jl*n;
a=il.'*h; 
b=a+h;
c=jl.'*h;
d=c+h;

q=zeros(1,ne);
L=@(u,v,x1,y1,x2,y2) lvsfun(x1+(x2-x1).*u+(y2-y1).*v,y1+(y2-y1).*u-(x2-x1).*v); % Change to new varibles (u,v)
Lu=@(u,v,x1,y1,x2,y2) (x2-x1).*lvsfunx(x1+(x2-x1).*u+(y2-y1).*v,y1+(y2-y1).*u-(x2-x1).*v)...
         +(y2-y1).*lvsfuny(x1+(x2-x1).*u+(y2-y1).*v,y1+(y2-y1).*u-(x2-x1).*v);
Lv=@(u,v,x1,y1,x2,y2) (y2-y1).*lvsfunx(x1+(x2-x1).*u+(y2-y1).*v,y1+(y2-y1).*u-(x2-x1).*v)...
         -(x2-x1).*lvsfuny(x1+(x2-x1).*u+(y2-y1).*v,y1+(y2-y1).*u-(x2-x1).*v);
vif=@(u,x1,y1,x2,y2) fzerov(@(u,v) L(u,v,x1,y1,x2,y2),@(u,v) Lv(u,v,x1,y1,x2,y2),u);
G=@(u,x1,y1,x2,y2) g(x1+(x2-x1).*u+(y2-y1).*vif(u,x1,y1,x2,y2),y1+(y2-y1).*u-(x2-x1).*vif(u,x1,y1,x2,y2))...
       .*sqrt((x2-x1).^2+(y2-y1).^2).*sqrt((Lu(u,vif(u,x1,y1,x2,y2),x1,y1,x2,y2)./Lv(u,vif(u,x1,y1,x2,y2),x1,y1,x2,y2)).^2+1);

% Type-1 or -2 elements: At the lower or upper part of the interface
ii=find(eli(2,:)==1 | eli(2,:)==2);
cv=upfemcv(@(u) vif(u,a(ii),eli(5,ii)',b(ii),eli(6,ii)'),zeros(length(ii),1),ones(length(ii),1));
q(ii)=quad1(@(u,l) G(u,a(ii(l)),eli(5,ii(l))',b(ii(l)),eli(6,ii(l))'),zeros(length(ii),1),ones(length(ii),1),p,cv);

% Type-3 or -4 elements: At the left or right part of the interface
ii=find(eli(2,:)==3 | eli(2,:)==4);
cv=upfemcv(@(u) vif(u,eli(3,ii)',c(ii),eli(4,ii)',d(ii)),zeros(length(ii),1),ones(length(ii),1));
q(ii)=quad1(@(u,l) G(u,eli(3,ii(l))',c(ii(l)),eli(4,ii(l))',d(ii(l))),zeros(length(ii),1),ones(length(ii),1),p,cv);

% Type-5 or 12 elements: At the lower-left corner, one vertex in Omega_1
%                        or the upper-right corner, three vertices in Omega_1
ii=find(eli(2,:)==5 | eli(2,:)==12);
cv=upfemcv(@(u) vif(u,eli(4,ii)',d(ii),b(ii),eli(6,ii)'),zeros(length(ii),1),ones(length(ii),1));
q(ii)=quad1(@(u,l) G(u,eli(4,ii(l))',d(ii(l)),b(ii(l)),eli(6,ii(l))'),zeros(length(ii),1),ones(length(ii),1),p,cv);

% Type-6 or 11 elements: At the lower-right corner, one vertex in Omega_1, 
%                        or the upper-left corner, three vertices in Omega_1.
ii=find(eli(2,:)==6 | eli(2,:)==11);
cv=upfemcv(@(u) vif(u,a(ii),eli(5,ii)',eli(4,ii)',d(ii)),zeros(length(ii),1),ones(length(ii),1));
q(ii)=quad1(@(u,l) G(u,a(ii(l)),eli(5,ii(l))',eli(4,ii(l))',d(ii(l))),zeros(length(ii),1),ones(length(ii),1),p,cv);

% Type-7 or 10 elements: At the upper-left corner, one vertex in Omega_1, 
%                        or the lower-right corner, three vertices in Omega_1
ii=find(eli(2,:)==7 | eli(2,:)==10);
cv=upfemcv(@(u) vif(u,eli(3,ii)',c(ii),b(ii),eli(6,ii)'),zeros(length(ii),1),ones(length(ii),1));
q(ii)=quad1(@(u,l) G(u,eli(3,ii(l))',c(ii(l)),b(ii(l)),eli(6,ii(l))'),zeros(length(ii),1),ones(length(ii),1),p,cv);

% Type-8 or 9 elements: At the upper-right corner, one vertex in Omega_1
%                       or the lower-left corner, three vertices in Omega_1
ii=find(eli(2,:)==8 | eli(2,:)==9);
cv=upfemcv(@(u) vif(u,a(ii),eli(5,ii)',eli(3,ii)',c(ii)),zeros(length(ii),1),ones(length(ii),1));
q(ii)=quad1(@(u,l) G(u,a(ii(l)),eli(5,ii(l))',eli(3,ii(l))',c(ii(l))),zeros(length(ii),1),ones(length(ii),1),p,cv);

% Type-13 or 21 elements
ii=find(eli(2,:)==13 | eli(2,:)==21);
cv=upfemcv(@(u) vif(u,eli(3,ii)',d(ii),eli(4,ii)',d(ii)),zeros(length(ii),1),ones(length(ii),1));
q(ii)=quad1(@(u,l) G(u,eli(3,ii(l))',d(ii(l)),eli(4,ii(l))',d(ii(l))),zeros(length(ii),1),ones(length(ii),1),p,cv);

% Type-14 or 22 elements
ii=find(eli(2,:)==14 | eli(2,:)==22);
cv1=upfemcv(@(u) vif(u,a(ii),eli(5,ii)',eli(3,ii)',c(ii)),zeros(length(ii),1),ones(length(ii),1));
cv2=upfemcv(@(u) vif(u,eli(4,ii)',c(ii),b(ii),eli(6,ii)'),zeros(length(ii),1),ones(length(ii),1));
q(ii)= quad1(@(u,l) G(u,a(ii(l)),eli(5,ii(l))',eli(3,ii(l))',c(ii(l))),zeros(length(ii),1),ones(length(ii),1),p,cv1)...
      +quad1(@(u,l) G(u,eli(4,ii(l))',c(ii(l)),b(ii(l)),eli(6,ii(l))'),zeros(length(ii),1),ones(length(ii),1),p,cv2);

% Type-15 or 23 elements
ii=find(eli(2,:)==15 | eli(2,:)==23);
cv=upfemcv(@(u) vif(u,eli(3,ii)',c(ii),eli(4,ii)',c(ii)),zeros(length(ii),1),ones(length(ii),1));
q(ii)=quad1(@(u,l) G(u,eli(3,ii(l))',c(ii(l)),eli(4,ii(l))',c(ii(l))),zeros(length(ii),1),ones(length(ii),1),p,cv);

% Type-16 or 24 elements
ii=find(eli(2,:)==16 | eli(2,:)==24);
cv1=upfemcv(@(u) vif(u,a(ii),eli(5,ii)',eli(3,ii)',d(ii)),zeros(length(ii),1),ones(length(ii),1));
cv2=upfemcv(@(u) vif(u,eli(4,ii)',d(ii),b(ii),eli(6,ii)'),zeros(length(ii),1),ones(length(ii),1));
q(ii)= quad1(@(u,l) G(u,a(ii(l)),eli(5,ii(l))',eli(3,ii(l))',d(ii(l))),zeros(length(ii),1),ones(length(ii),1),p,cv1)...
      +quad1(@(u,l) G(u,eli(4,ii(l))',d(ii(l)),b(ii(l)),eli(6,ii(l))'),zeros(length(ii),1),ones(length(ii),1),p,cv2);

% Type-17 or 25 elements
ii=find(eli(2,:)==17 | eli(2,:)==25);
cv=upfemcv(@(u) vif(u,b(ii),eli(5,ii)',b(ii),eli(6,ii)'),zeros(length(ii),1),ones(length(ii),1));
q(ii)=quad1(@(u,l) G(u,b(ii(l)),eli(5,ii(l))',b(ii(l)),eli(6,ii(l))'),zeros(length(ii),1),ones(length(ii),1),p,cv);

% Type-18 or 26 elements
ii=find(eli(2,:)==18 | eli(2,:)==26);
cv1=upfemcv(@(u) vif(u,eli(3,ii)',c(ii),a(ii),eli(5,ii)'),zeros(length(ii),1),ones(length(ii),1));
cv2=upfemcv(@(u) vif(u,a(ii),eli(6,ii)',eli(4,ii)',d(ii)),zeros(length(ii),1),ones(length(ii),1));
q(ii)= quad1(@(u,l) G(u,eli(3,ii(l))',c(ii(l)),a(ii(l)),eli(5,ii(l))'),zeros(length(ii),1),ones(length(ii),1),p,cv1)...
      +quad1(@(u,l) G(u,a(ii(l)),eli(6,ii(l))',eli(4,ii(l))',d(ii(l))),zeros(length(ii),1),ones(length(ii),1),p,cv2);

% Type-19 or 27 elements
ii=find(eli(2,:)==19 | eli(2,:)==27);
cv=upfemcv(@(u) vif(u,a(ii),eli(5,ii)',a(ii),eli(6,ii)'),zeros(length(ii),1),ones(length(ii),1));
q(ii)=quad1(@(u,l) G(u,a(ii(l)),eli(5,ii(l))',a(ii(l)),eli(6,ii(l))'),zeros(length(ii),1),ones(length(ii),1),p,cv);

% Type-20 or 28 elements
ii=find(eli(2,:)==20 | eli(2,:)==28);
cv1=upfemcv(@(u) vif(u,eli(3,ii)',c(ii),b(ii),eli(5,ii)'),zeros(length(ii),1),ones(length(ii),1));
cv2=upfemcv(@(u) vif(u,b(ii),eli(6,ii)',eli(4,ii)',d(ii)),zeros(length(ii),1),ones(length(ii),1));
q(ii)= quad1(@(u,l) G(u,eli(3,ii(l))',c(ii(l)),b(ii(l)),eli(5,ii(l))'),zeros(length(ii),1),ones(length(ii),1),p,cv1)...
      +quad1(@(u,l) G(u,b(ii(l)),eli(6,ii(l))',eli(4,ii(l))',d(ii(l))),zeros(length(ii),1),ones(length(ii),1),p,cv2);
