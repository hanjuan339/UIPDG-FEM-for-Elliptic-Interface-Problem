function [err1,err0,errinf,errflux,u,el1,el2,eli,ltg,cn,en1_2,en1_4,en2_2,en2_4]=ex1fun(a,theta,n,p,beta,gamma0,gamma1)
%Errors in H^1 semi-norm, L^2 norm, and L^infty norm for Exmaple 1.
%Exact solution and their partial derivatives in Omega_1 and Oemga_2
ue1=@(x,y) 1/a(1)*exp(x.*y);
ue1x=@(x,y) 1/a(1)*y.*exp(x.*y);
ue1y=@(x,y) 1/a(1)*x.*exp(x.*y);
ue2=@(x,y) 1/a(2)*sin(pi*x).*sin(pi*y);
ue2x=@(x,y) 1/a(2)*pi*cos(pi*x).*sin(pi*y);
ue2y=@(x,y) 1/a(2)*pi*sin(pi*x).*cos(pi*y);



%Leveset function and its partive derivatives.


% r0=1/10; b0=sqrt(10); %r=r0*(b0+cos(6*th)) % flower-6
%                       %Known issue with r0=0.1;b0=3: the line y=0.7 tangent the interface and the intersect point (0.5,0.7) can not be approximated well by fzerox
% %r0=sqrt(2)/14; b0=3;
% L=@(x,y) ((x-0.5).^2+(y-0.5).^2-0.01).^2.*((x-0.5).^2+(y-0.5).^2<0.01);
% Lx=@(x,y) 4*(x-0.5).*((x-0.5).^2+(y-0.5).^2-0.01).*((x-0.5).^2+(y-0.5).^2<0.01);
% Ly=@(x,y) 4*(y-0.5).*((x-0.5).^2+(y-0.5).^2-0.01).*((x-0.5).^2+(y-0.5).^2<0.01);
% lvsfun=@(x,y) ((x-0.5).^2+(y-0.5).^2).^(7/2)-r0*(b0*((x-0.5).^2+(y-0.5).^2).^3+(x-0.5).^6-(y-0.5).^6-15*(x-0.5).^4.*(y-0.5).^2+15*(x-0.5).^2.*(y-0.5).^4)-L(x,y); 
% lvsfunx=@(x,y) (x-0.5).*(7*((x-0.5).^2+(y-0.5).^2).^(5/2)-6*r0*(b0*((x-0.5).^2+(y-0.5).^2).^2+(x-0.5).^4-10*(x-0.5).^2.*(y-0.5).^2+5*(y-0.5).^4))-Lx(x,y);
% lvsfuny=@(x,y) (y-0.5).*(7*((x-0.5).^2+(y-0.5).^2).^(5/2)-6*r0*(b0*((x-0.5).^2+(y-0.5).^2).^2-(y-0.5).^4+10*(x-0.5).^2.*(y-0.5).^2-5*(x-0.5).^4))-Ly(x,y);



lvsfun=@(x,y) ((2*x-1).^2+(2*y-1).^2).^3-1/2*sqrt((2*x-1).^2+(2*y-1).^2).^5-(5*(2*x-1).^4.*(2*y-1)-10*(2*x-1).^2.*(2*y-1).^3+(2*y-1).^5)/7; %flower-5
lvsfunx=@(x,y) 5*(1-2*x).*sqrt((2*x-1).^2+(2*y-1).^2).^3-12*(1-2*x).*((1-2*x).^2+(1-2*y).^2).^2+(40*(1-2*x).^3.*(2*y-1)-40*(1-2*x).*(2*y-1).^3)/7;
lvsfuny=@(x,y) 5*(1-2*y).*sqrt((1-2*x).^2+(1-2*y).^2).^3-12*(1-2*y).*((1-2*x).^2+(1-2*y).^2).^2+(-10*(1-2*x).^4+60*(1-2*x).^2.*(1-2*y).^2-10*(1-2*y).^4)/7;

% lvsfun=@(x,y) (2*((2*x-1/2).^2+(2*y-1).^2)-2*x+1/2).^2-((2*x-1/2).^2+(2*y-1).^2)+1/10;%Kidney-shaped interface
% lvsfunx=@(x,y) 4*(8*x-3).*(2*((2*x-1/2).^2+(2*y-1).^2)-2*x+1/2)-2*(4*x-1);
% lvsfuny=@(x,y) 16*(2*y-1).*(2*((2*x-1/2).^2+(2*y-1).^2)-2*x+1/2)-4*(2*y-1);

% lvsfun=@(x,y) ((2*x-1).^2+(2*y-1).^2).^2-1/2*sqrt((2*x-1).^2+(2*y-1).^2).^3-(3*(2*x-1).^2.*(2*y-1)-(2*y-1).^3)/7;%flower-3
% lvsfunx=@(x,y) 8*(2*x-1).*((2*x-1).^2+(2*y-1).^2)-3*(2*x-1).*sqrt((2*x-1).^2+(2*y-1).^2)-12*(2*x-1).*(2*y-1)/7;
% lvsfuny=@(x,y) 8*(2*y-1).*((2*x-1).^2+(2*y-1).^2)-3*(2*y-1).*sqrt((2*x-1).^2+(2*y-1).^2)-6*((2*x-1).^2-(2*y-1).^2)/7;

% lvsfun=@(x,y) ((2*x-1).^2+(2*y-1).^2).^4-1/2*sqrt((2*x-1).^2+(2*y-1).^2).^7-(7*(2*x-1).^6.*(2*y-1)+21*(2*x-1).^2.*(2*y-1).^5-35*(2*x-1).^4.*(2*y-1).^3-(2*y-1).^7)/7;%flower-7
% lvsfunx=@(x,y) 16*(2*x-1).*((2*x-1).^2+(2*y-1).^2).^3-7*(2*x-1).*sqrt((2*x-1).^2+(2*y-1).^2).^5-(12*(2*x-1).^5.*(2*y-1)+12*(2*x-1).*(2*y-1).^5-40*(2*x-1).^3.*(2*y-1).^3);
% lvsfuny=@(x,y) 16*(2*y-1).*((2*x-1).^2+(2*y-1).^2).^3-7*(2*y-1).*sqrt((2*x-1).^2+(2*y-1).^2).^5-(2*(2*x-1).^6+30*(2*x-1).^2.*(2*y-1).^4-30*(2*x-1).^4.*(2*y-1).^2-2*(2*y-1).^6);

% lvsfun=@(x,y) (x-0.5).^2+(y-0.5).^2-1/16; %cicle
% lvsfunx=@(x,y) (2*x-1);
% lvsfuny=@(x,y) (2*y-1);


%Data
f1=@(x,y) -(x.^2+y.^2).*exp(x.*y)+theta(1)/a(1)*exp(x.*y);
f2=@(x,y) 2*pi^2*sin(pi*x).*sin(pi*y)+theta(2)/a(2)*sin(pi*x).*sin(pi*y);

gd=@(x,y) ue1(x,y)-ue2(x,y);
gn=@(x,y) ((a(1)*ue1x(x,y)-a(2)*ue2x(x,y)).*lvsfunx(x,y)+(a(1)*ue1y(x,y)-a(2)*ue2y(x,y)).*lvsfuny(x,y))./sqrt(lvsfunx(x,y).^2+lvsfuny(x,y).^2);

%Solution
[u,el1,el2,eli,ltg,cn,en1_2,en1_4,en2_2,en2_4,~,~]=uipdgfem(a,theta,f1,f2,gd,gn,lvsfun,lvsfunx,lvsfuny,n,p,beta,gamma0,gamma1);

%Errors
[err1,err0,errinf,errflux]=upfemerr(u,ue1,ue2,f1,f2,lvsfun,el1,el2,eli,ltg,n,p,lvsfunx,lvsfuny,ue1x,ue1y,ue2x,ue2y,a);
