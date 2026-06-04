function q=upfemitfquad2d(g,lvsfun,lvsfunx,lvsfuny,eli,n,p,sgn)
%Integration over each of the two parts of square eleemnts seperated by a interface.
%   g: integrand
%   lvsfunx: the partial derivative of the lvsfun w.r.t. x
%   lvsfuny: the partial derivative of the lvsfun w.r.t. y
%   eli: Interface elements.
%   q: the integrations of g over the parts in Omega_1 if sgn=-1, or in Omega_2 if sgn=1.

ne=size(eli,2); % Number of interface elements
h=1/n;
G=@(x,y) g(y,x);

% Conners of interface elements
el=eli(1,:);
jl=floor((el-1)/n);
il=el-1-jl*n;
a=il.'*h; 
b=a+h;
c=jl.'*h;
d=c+h;

q=zeros(1,ne);
% Type-1 elements: The interface cuts the left and right edges, the upper edge locates in Omega_1
ii=find(eli(2,:)==1);
yif=@(x,l) fzeroy(lvsfun,lvsfuny,x,(eli(5,ii(l))+eli(6,ii(l)))'/2);
if sgn<0
    q(ii)=quad2(g,a(ii),b(ii),yif,d(ii),p);
else
    q(ii)=quad2(g,a(ii),b(ii),c(ii),yif,p);
end
% Type-2 elements: The interface cuts the left and right edges, the lower edge in Omega_1
ii=find(eli(2,:)==2);
yif=@(x,l) fzeroy(lvsfun,lvsfuny,x,(eli(5,ii(l))+eli(6,ii(l)))'/2);
if sgn<0
    q(ii)=quad2(g,a(ii),b(ii),c(ii),yif,p);
else
    q(ii)=quad2(g,a(ii),b(ii),yif,d(ii),p);
end
% Type-3 elements: The interface cuts the lower and upper edges, the right edge in Omega_1
%                  Evaluate integrals by interchanging x and y
ii=find(eli(2,:)==3);
xif=@(y,l) fzerox(lvsfun,lvsfunx,(eli(3,ii(l))+eli(4,ii(l)))'/2,y);
if sgn<0
    q(ii)=quad2(G,c(ii),d(ii),xif,b(ii),p);
else
    q(ii)=quad2(G,c(ii),d(ii),a(ii),xif,p);
end
% Type-4 elements: The interface cuts the lower and upper edges, the left edge in Omega_1
%                  Evaluate integrals by interchanging x and y
ii=find(eli(2,:)==4);
xif=@(y,l) fzerox(lvsfun,lvsfunx,(eli(3,ii(l))+eli(4,ii(l)))'/2,y);
if sgn<0
    q(ii)=quad2(G,c(ii),d(ii),a(ii),xif,p);
else
    q(ii)=quad2(G,c(ii),d(ii),xif,b(ii),p);
end
% Type-5 elements: The interface cuts the upper and right edges, the upper-right conner in Omega_1
ii=find(eli(2,:)==5);
% jj=atan(lvsfuny(eli(4,ii)',d(ii))./lvsfunx(eli(4,ii)',d(ii)))<max(pi/10,atan(lvsfunx(b(ii),eli(6,ii)')./lvsfuny(b(ii),eli(6,ii)')))...
%    | abs(atan(lvsfuny(b(ii),eli(6,ii)')./lvsfunx(b(ii),eli(6,ii)')))<pi/10;
jj=atan(lvsfuny(eli(4,ii)',d(ii))./lvsfunx(eli(4,ii)',d(ii)))<atan(lvsfunx(b(ii),eli(6,ii)')./lvsfuny(b(ii),eli(6,ii)'));
i1=ii(~jj);
i2=ii(jj);
xc=eli(4,i1)'; %fzerox(lvsfun,lvsfunx,a(i1)+0.5*h,d(i1));
yc=eli(6,i2)'; %fzeroy(lvsfun,lvsfuny,b(i2),c(i2)+0.5*h);
xif=@(y,l) fzerox(lvsfun,lvsfunx,(eli(4,i2(l))'+b(i2(l)))/2,y);
yif=@(x,l) fzeroy(lvsfun,lvsfuny,x,(eli(6,i1(l))'+d(i1(l)))/2);
if sgn<0
    q(i1)=quad2(g,xc,b(i1),yif,d(i1),p);
    q(i2)=quad2(G,yc,d(i2),xif,b(i2),p);
else
    q(i1)=quad2(g,a(i1),xc,c(i1),d(i1),p)+quad2(g,xc,b(i1),c(i1),yif,p);
    q(i2)=quad2(G,c(i2),yc,a(i2),b(i2),p)+quad2(G,yc,d(i2),a(i2),xif,p);
end
% Type-6 elements: The interface cuts the upper and left edges, the upper-left conner in Omega_1
ii=find(eli(2,:)==6);
% jj=atan(lvsfuny(eli(4,ii)',d(ii))./lvsfunx(eli(4,ii)',d(ii)))>min(-pi/10,atan(lvsfunx(a(ii),eli(5,ii)')./lvsfuny(a(ii),eli(5,ii)')))...
%    | abs(atan(lvsfuny(a(ii),eli(5,ii)')./lvsfunx(a(ii),eli(5,ii)')))<pi/10;
jj=atan(-lvsfuny(eli(4,ii)',d(ii))./lvsfunx(eli(4,ii)',d(ii)))<atan(-lvsfunx(a(ii),eli(5,ii)')./lvsfuny(a(ii),eli(5,ii)'));
i1=ii(~jj);
i2=ii(jj);
xc=eli(4,i1)'; %fzerox(lvsfun,lvsfunx,a(i1)+0.5*h,d(i1));
yc=eli(5,i2)'; %fzeroy(lvsfun,lvsfuny,a(i2),c(i2)+0.5*h);
xif=@(y,l) fzerox(lvsfun,lvsfunx,(eli(4,i2(l))'+a(i2(l)))/2,y);
yif=@(x,l) fzeroy(lvsfun,lvsfuny,x,(eli(5,i1(l))'+d(i1(l)))/2);
if sgn<0
    q(i1)=quad2(g,a(i1),xc,yif,d(i1),p);
    q(i2)=quad2(G,yc,d(i2),a(i2),xif,p);
else
    q(i1)=quad2(g,xc,b(i1),c(i1),d(i1),p)+quad2(g,a(i1),xc,c(i1),yif,p);
    q(i2)=quad2(G,c(i2),yc,a(i2),b(i2),p)+quad2(G,yc,d(i2),xif,b(i2),p);
end
% Type-7 elements: The interface cuts the lower and right edges, the lower-right conner in Omega_1
ii=find(eli(2,:)==7);
% jj=atan(lvsfuny(eli(3,ii)',c(ii))./lvsfunx(eli(3,ii)',c(ii)))>min(-pi/10,atan(lvsfunx(b(ii),eli(6,ii)')./lvsfuny(b(ii),eli(6,ii)')))...
%    | abs(atan(lvsfuny(b(ii),eli(6,ii)')./lvsfunx(b(ii),eli(6,ii)')))<pi/10;
jj=atan(-lvsfuny(eli(3,ii)',c(ii))./lvsfunx(eli(3,ii)',c(ii)))<atan(-lvsfunx(b(ii),eli(6,ii)')./lvsfuny(b(ii),eli(6,ii)'));
i1=ii(~jj);
i2=ii(jj);
xc=eli(3,i1)'; %fzerox(lvsfun,lvsfunx,a(i1)+0.5*h,c(i1));
yc=eli(6,i2)'; %fzeroy(lvsfun,lvsfuny,b(i2),c(i2)+0.5*h);
xif=@(y,l) fzerox(lvsfun,lvsfunx,(eli(3,i2(l))'+b(i2(l)))/2,y);
yif=@(x,l) fzeroy(lvsfun,lvsfuny,x,(eli(6,i1(l))'+c(i1(l)))/2);
if sgn<0
    q(i1)=quad2(g,xc,b(i1),c(i1),yif,p);
    q(i2)=quad2(G,c(i2),yc,xif,b(i2),p);
else
    q(i1)=quad2(g,a(i1),xc,c(i1),d(i1),p)+quad2(g,xc,b(i1),yif,d(i1),p);
    q(i2)=quad2(G,yc,d(i2),a(i2),b(i2),p)+quad2(G,c(i2),yc,a(i2),xif,p);
end
% Type-8 elements: The interface cuts the lower and left edges, the lower-left conner in Omega_1
ii=find(eli(2,:)==8);
% jj=atan(lvsfuny(eli(3,ii)',c(ii))./lvsfunx(eli(3,ii)',c(ii)))<max(pi/10,atan(lvsfunx(a(ii),eli(5,ii)')./lvsfuny(a(ii),eli(5,ii)')))...
%    | abs(atan(lvsfuny(a(ii),eli(5,ii)')./lvsfunx(a(ii),eli(5,ii)')))<pi/10;
jj=atan(lvsfuny(eli(3,ii)',c(ii))./lvsfunx(eli(3,ii)',c(ii)))<atan(lvsfunx(a(ii),eli(5,ii)')./lvsfuny(a(ii),eli(5,ii)'));
i1=ii(~jj);
i2=ii(jj);
xc=eli(3,i1)'; %fzerox(lvsfun,lvsfunx,a(i1)+0.5*h,c(i1));
yc=eli(5,i2)'; %fzeroy(lvsfun,lvsfuny,a(i2),c(i2)+0.5*h);
xif=@(y,l) fzerox(lvsfun,lvsfunx,(eli(3,i2(l))'+a(i2(l)))/2,y);
yif=@(x,l) fzeroy(lvsfun,lvsfuny,x,(eli(5,i1(l))'+c(i1(l)))/2);
if sgn<0
    q(i1)=quad2(g,a(i1),xc,c(i1),yif,p);
    q(i2)=quad2(G,c(i2),yc,a(i2),xif,p);
else
    q(i1)=quad2(g,xc,b(i1),c(i1),d(i1),p)+quad2(g,a(i1),xc,yif,d(i1),p);
    q(i2)=quad2(G,yc,d(i2),a(i2),b(i2),p)+quad2(G,c(i2),yc,xif,b(i2),p);
end
% Type-9 elements: The interface cuts the lower and left edges, the lower-left conner in Omega_2
ii=find(eli(2,:)==9);
% jj=atan(lvsfuny(eli(3,ii)',c(ii))./lvsfunx(eli(3,ii)',c(ii)))<max(pi/10,atan(lvsfunx(a(ii),eli(5,ii)')./lvsfuny(a(ii),eli(5,ii)')))...
%    | abs(atan(lvsfuny(a(ii),eli(5,ii)')./lvsfunx(a(ii),eli(5,ii)')))<pi/10;
jj=atan(lvsfuny(eli(3,ii)',c(ii))./lvsfunx(eli(3,ii)',c(ii)))<atan(lvsfunx(a(ii),eli(5,ii)')./lvsfuny(a(ii),eli(5,ii)'));
i1=ii(~jj);
i2=ii(jj);
xc=eli(3,i1)'; %fzerox(lvsfun,lvsfunx,a(i1)+0.5*h,c(i1));
yc=eli(5,i2)'; %fzeroy(lvsfun,lvsfuny,a(i2),c(i2)+0.5*h);
xif=@(y,l) fzerox(lvsfun,lvsfunx,(eli(3,i2(l))'+a(i2(l)))/2,y);
yif=@(x,l) fzeroy(lvsfun,lvsfuny,x,(eli(5,i1(l))'+c(i1(l)))/2);
if sgn<0
    q(i1)=quad2(g,xc,b(i1),c(i1),d(i1),p)+quad2(g,a(i1),xc,yif,d(i1),p);
    q(i2)=quad2(G,yc,d(i2),a(i2),b(i2),p)+quad2(G,c(i2),yc,xif,b(i2),p);
else
    q(i1)=quad2(g,a(i1),xc,c(i1),yif,p);
    q(i2)=quad2(G,c(i2),yc,a(i2),xif,p);
end
% Type-10 elements: The interface cuts the lower and right edges, the lower-right conner in Omega_2
ii=find(eli(2,:)==10);
% jj=atan(lvsfuny(eli(3,ii)',c(ii))./lvsfunx(eli(3,ii)',c(ii)))>min(-pi/10,atan(lvsfunx(b(ii),eli(6,ii)')./lvsfuny(b(ii),eli(6,ii)')))...
%    | abs(atan(lvsfuny(b(ii),eli(6,ii)')./lvsfunx(b(ii),eli(6,ii)')))<pi/10;
jj=atan(-lvsfuny(eli(3,ii)',c(ii))./lvsfunx(eli(3,ii)',c(ii)))<atan(-lvsfunx(b(ii),eli(6,ii)')./lvsfuny(b(ii),eli(6,ii)'));
i1=ii(~jj);
i2=ii(jj);
xc=eli(3,i1)'; %fzerox(lvsfun,lvsfunx,a(i1)+0.5*h,c(i1));
yc=eli(6,i2)'; %fzeroy(lvsfun,lvsfuny,b(i2),c(i2)+0.5*h);
xif=@(y,l) fzerox(lvsfun,lvsfunx,(eli(3,i2(l))'+b(i2(l)))/2,y);
yif=@(x,l) fzeroy(lvsfun,lvsfuny,x,(eli(6,i1(l))'+c(i1(l)))/2);
if sgn<0
    q(i1)=quad2(g,a(i1),xc,c(i1),d(i1),p)+quad2(g,xc,b(i1),yif,d(i1),p);
    q(i2)=quad2(G,yc,d(i2),a(i2),b(i2),p)+quad2(G,c(i2),yc,a(i2),xif,p);
else
    q(i1)=quad2(g,xc,b(i1),c(i1),yif,p);
    q(i2)=quad2(G,c(i2),yc,xif,b(i2),p);
end
% Type-11 elements: The interface cuts the upper and left edges, the upper-left conner in Omega_2
ii=find(eli(2,:)==11);
% jj=atan(lvsfuny(eli(4,ii)',d(ii))./lvsfunx(eli(4,ii)',d(ii)))>min(-pi/10,atan(lvsfunx(a(ii),eli(5,ii)')./lvsfuny(a(ii),eli(5,ii)')))...
%    | abs(atan(lvsfuny(a(ii),eli(5,ii)')./lvsfunx(a(ii),eli(5,ii)')))<pi/10;
jj=atan(-lvsfuny(eli(4,ii)',d(ii))./lvsfunx(eli(4,ii)',d(ii)))<atan(-lvsfunx(a(ii),eli(5,ii)')./lvsfuny(a(ii),eli(5,ii)'));
i1=ii(~jj);
i2=ii(jj);
xc=eli(4,i1)'; %fzerox(lvsfun,lvsfunx,a(i1)+0.5*h,d(i1));
yc=eli(5,i2)'; %fzeroy(lvsfun,lvsfuny,a(i2),c(i2)+0.5*h);
xif=@(y,l) fzerox(lvsfun,lvsfunx,(eli(4,i2(l))'+a(i2(l)))/2,y);
yif=@(x,l) fzeroy(lvsfun,lvsfuny,x,(eli(5,i1(l))'+d(i1(l)))/2);
if sgn<0
    q(i1)=quad2(g,xc,b(i1),c(i1),d(i1),p)+quad2(g,a(i1),xc,c(i1),yif,p);
    q(i2)=quad2(G,c(i2),yc,a(i2),b(i2),p)+quad2(G,yc,d(i2),xif,b(i2),p);
else
    q(i1)=quad2(g,a(i1),xc,yif,d(i1),p);
    q(i2)=quad2(G,yc,d(i2),a(i2),xif,p);
end
% Type-12 elements: The interface cuts the upper and right edges, the upper-right conner in Omega_2
ii=find(eli(2,:)==12);
% jj=atan(lvsfuny(eli(4,ii)',d(ii))./lvsfunx(eli(4,ii)',d(ii)))<max(pi/10,atan(lvsfunx(b(ii),eli(6,ii)')./lvsfuny(b(ii),eli(6,ii)')))...
%    | abs(atan(lvsfuny(b(ii),eli(6,ii)')./lvsfunx(b(ii),eli(6,ii)')))<pi/10;
jj=atan(lvsfuny(eli(4,ii)',d(ii))./lvsfunx(eli(4,ii)',d(ii)))<atan(lvsfunx(b(ii),eli(6,ii)')./lvsfuny(b(ii),eli(6,ii)'));
i1=ii(~jj);
i2=ii(jj);
xc=eli(4,i1)'; %fzerox(lvsfun,lvsfunx,a(i1)+0.5*h,d(i1));
yc=eli(6,i2)'; %fzeroy(lvsfun,lvsfuny,b(i2),c(i2)+0.5*h);
xif=@(y,l) fzerox(lvsfun,lvsfunx,(eli(4,i2(l))'+b(i2(l)))/2,y);
yif=@(x,l) fzeroy(lvsfun,lvsfuny,x,(eli(6,i1(l))'+d(i1(l)))/2);
if sgn<0
    q(i1)=quad2(g,a(i1),xc,c(i1),d(i1),p)+quad2(g,xc,b(i1),c(i1),yif,p);
    q(i2)=quad2(G,c(i2),yc,a(i2),b(i2),p)+quad2(G,yc,d(i2),a(i2),xif,p);
else
    q(i1)=quad2(g,xc,b(i1),yif,d(i1),p);
    q(i2)=quad2(G,yc,d(i2),xif,b(i2),p);
end

% Type-13 or 21 elements: The interface cuts the upper edge twice, four conners in Omega_2 (13) or Omega_1 (21)
ii=find((sgn<0 & eli(2,:)==13) | (sgn>0 & eli(2,:)==21));
xl=eli(3,ii)';
xr=eli(4,ii)';
yif=@(x,l) fzeroy(lvsfun,lvsfuny,x,d(ii(l)));
q(ii)=quad2(g,xl,xr,yif,d(ii),p);
ii=find((sgn>0 & eli(2,:)==13) | (sgn<0 & eli(2,:)==21));
xl=eli(3,ii)';
xr=eli(4,ii)';
yif=@(x,l) fzeroy(lvsfun,lvsfuny,x,d(ii(l)));
q(ii)=quad2(g,a(ii),xl,c(ii),d(ii),p)+quad2(g,xl,xr,c(ii),yif,p)+quad2(g,xr,b(ii),c(ii),d(ii),p);

% Type-14 or 22 elements: The interface cuts the lower edge twice and the left and right edges once, respectively, the upper edge in Omega_1 (14) or Omega_2 (22)
ii=find((sgn<0 & eli(2,:)==14) | (sgn>0 & eli(2,:)==22));
xl=eli(3,ii)';
xr=eli(4,ii)';
yif=@(x,l) fzeroy(lvsfun,lvsfuny,x,c(ii(l)));
q(ii)=quad2(g,a(ii),xl,yif,d(ii),p)+quad2(g,xl,xr,c(ii),d(ii),p)+quad2(g,xr,b(ii),yif,d(ii),p);
ii=find((sgn>0 & eli(2,:)==14) | (sgn<0 & eli(2,:)==22));
xl=eli(3,ii)';
xr=eli(4,ii)';
yif=@(x,l) fzeroy(lvsfun,lvsfuny,x,c(ii(l)));
q(ii)=quad2(g,a(ii),xl,c(ii),yif,p)+quad2(g,xr,b(ii),c(ii),yif,p);

% Type-15 or 23 elements: The interface cuts the lower edge twice, four conners in Omega_2 (15) or Omega_1 (23)
ii=find((sgn<0 & eli(2,:)==15) | (sgn>0 & eli(2,:)==23));
xl=eli(3,ii)';
xr=eli(4,ii)';
yif=@(x,l) fzeroy(lvsfun,lvsfuny,x,c(ii(l)));
q(ii)=quad2(g,xl,xr,c(ii),yif,p);
ii=find((sgn>0 & eli(2,:)==15) | (sgn<0 & eli(2,:)==23));
xl=eli(3,ii)';
xr=eli(4,ii)';
yif=@(x,l) fzeroy(lvsfun,lvsfuny,x,c(ii(l)));
q(ii)=quad2(g,a(ii),xl,c(ii),d(ii),p)+quad2(g,xl,xr,yif,d(ii),p)+quad2(g,xr,b(ii),c(ii),d(ii),p);

% Type-16 or 24 elements: The interface cuts the upper edge twice and the left and right edges once, respectively, the lower edge in Omega_1 (16) or Omega_2 (21)
ii=find((sgn<0 & eli(2,:)==16) | (sgn>0 & eli(2,:)==24));
xl=eli(3,ii)';
xr=eli(4,ii)';
yif=@(x,l) fzeroy(lvsfun,lvsfuny,x,d(ii(l)));
q(ii)=quad2(g,a(ii),xl,c(ii),yif,p)+quad2(g,xl,xr,c(ii),d(ii),p)+quad2(g,xr,b(ii),c(ii),yif,p);
ii=find((sgn>0 & eli(2,:)==16) | (sgn<0 & eli(2,:)==24));
xl=eli(3,ii)';
xr=eli(4,ii)';
yif=@(x,l) fzeroy(lvsfun,lvsfuny,x,d(ii(l)));
q(ii)=quad2(g,a(ii),xl,yif,d(ii),p)+quad2(g,xr,b(ii),yif,d(ii),p);

% Type-17 or 25 elements: The interface cuts the right edge twice, four conners in Omega_2 (17) or Omega_1 (25)
ii=find((sgn<0 & eli(2,:)==17) | (sgn>0 & eli(2,:)==25));
yl=eli(5,ii)';
yu=eli(6,ii)';
xif=@(y,l) fzerox(lvsfun,lvsfunx,b(ii(l)),y);
q(ii)=quad2(G,yl,yu,xif,b(ii),p);
ii=find((sgn>0 & eli(2,:)==17) | (sgn<0 & eli(2,:)==25));
yl=eli(5,ii)';
yu=eli(6,ii)';
xif=@(y,l) fzerox(lvsfun,lvsfunx,b(ii(l)),y);
q(ii)=quad2(G,c(ii),yl,a(ii),b(ii),p)+quad2(G,yl,yu,a(ii),xif,p)+quad2(G,yu,d(ii),a(ii),b(ii),p);

% Type-18 or 26 elements: The interface cuts the left edge twice and the lower and upper edges once, respectively, the right edge in Omega_1 (18) or Omega_2 (26)
ii=find((sgn<0 & eli(2,:)==18) | (sgn>0 & eli(2,:)==26));
yl=eli(5,ii)';
yu=eli(6,ii)';
xif=@(y,l) fzerox(lvsfun,lvsfunx,a(ii(l)),y);
q(ii)=quad2(G,c(ii),yl,xif,b(ii),p)+quad2(G,yl,yu,a(ii),b(ii),p)+quad2(G,yu,d(ii),xif,b(ii),p);
ii=find((sgn>0 & eli(2,:)==18) | (sgn<0 & eli(2,:)==26));
yl=eli(5,ii)';
yu=eli(6,ii)';
xif=@(y,l) fzerox(lvsfun,lvsfunx,a(ii(l)),y);
q(ii)=quad2(G,c(ii),yl,a(ii),xif,p)+quad2(G,yu,d(ii),a(ii),xif,p);

% Type-19 or 27 elements: The interface cuts the left edge twice, four conners in Omega_2 (19) or Omega_1 (27)
ii=find((sgn<0 & eli(2,:)==19) | (sgn>0 & eli(2,:)==27));
yl=eli(5,ii)';
yu=eli(6,ii)';
xif=@(y,l) fzerox(lvsfun,lvsfunx,a(ii(l)),y);
q(ii)=quad2(G,yl,yu,a(ii),xif,p);
ii=find((sgn>0 & eli(2,:)==19) | (sgn<0 & eli(2,:)==27));
yl=eli(5,ii)';
yu=eli(6,ii)';
xif=@(y,l) fzerox(lvsfun,lvsfunx,a(ii(l)),y);
q(ii)=quad2(G,c(ii),yl,a(ii),b(ii),p)+quad2(G,yl,yu,xif,b(ii),p)+quad2(G,yu,d(ii),a(ii),b(ii),p);

% Type-20 or 28 elements: The interface cuts the right edge twice and the lower and upper edges once, respectively, the left edge in Omega_1 (20) or Omega_2 (28)
ii=find((sgn<0 & eli(2,:)==20) | (sgn>0 & eli(2,:)==28));
yl=eli(5,ii)';
yu=eli(6,ii)';
xif=@(y,l) fzerox(lvsfun,lvsfunx,b(ii(l)),y);
q(ii)=quad2(G,c(ii),yl,a(ii),xif,p)+quad2(G,yl,yu,a(ii),b(ii),p)+quad2(G,yu,d(ii),a(ii),xif,p);
ii=find((sgn>0 & eli(2,:)==20) | (sgn<0 & eli(2,:)==28));
yl=eli(5,ii)';
yu=eli(6,ii)';
xif=@(y,l) fzerox(lvsfun,lvsfunx,b(ii(l)),y);
q(ii)=quad2(G,c(ii),yl,xif,b(ii),p)+quad2(G,yu,d(ii),xif,b(ii),p);
