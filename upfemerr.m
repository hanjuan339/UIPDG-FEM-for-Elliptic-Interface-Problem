function [err1,err0,errinf,errflux]=upfemerr(u,ue1,ue2,f1,f2,lvsfun,el1,el2,eli,ltg,n,p,lvsfunx,lvsfuny,ue1x,ue1y,ue2x,ue2y,a)
%Evaluate relative errors between the UPFE solution and the exact solution.
%
%   u: Nodal values of the discrete solution.
%   el1: elements inside the interface (in Omega_1 where the lvsfun<0)
%   el2: elements outside the interface (in Omega_2 where the lvsfun>0)
%   eli: interface elements.
%   ltgmapl,l=1,2: the local to global mapping in Omega_l. ltgmaps is 
%           a matrix with (p+1)^2 rows. The m-th column gives the DOFs 
%           in the ell(m)-th element, in precise, ltgmaps(j,m) is the index 
%           to the golbal DOF of the j-th DOF in the m-th element.
%   ltgmapil,l=1,2: the local to global mapping for elements crossed by 
%           the interface. Note that elements crossed by the interface are
%           assigned two set of DOFs related to Omega_l. l=1, 2,
%           respectively.
%   n: h=1/n is the mesh size
%   p: the degree of element polynomials (in each x and y, respectively).
%   normflag=1: Relative energy norm
%           =0: Relative L^2 norm
%           =inf: Relative L^infinity norm.
%   uel,l=1,2: exact solution in Omega_l.
%   uelx,uely, l=1,2: Gradients of uel. No needed for L^2 and L^infty norms.
%   err(l): the error in Omega_l.

ltgmap1=ltg{1};
ltgmap2=ltg{2};
ltgmapi1=ltg{3};
ltgmapi2=ltg{4};


h=1/n;
xn=linspace(0,1,p+1).'; %Nodal points, equal distances

% Relative error in L^infinity norm.
ui=zeros(size(u)); %Interpolant of the exact solution.
el1i=[el1 eli(1,:)];
j=floor((el1i-1)/n);
i=el1i-1-j*n;
x=[];
y=[];
for l=1:p+1
    for m=1:p+1
        x=[x;(xn(m)+i)*h];
        y=[y;(xn(l)+j)*h];
    end
end
uh1=u([ltgmap1 ltgmapi1]);
u1=ue1(x,y);
ui([ltgmap1 ltgmapi1])=u1;
uh1=uh1.*(lvsfun(x,y)<=eps);%Get rid of DOFs in Omega_2.
u1=u1.*(lvsfun(x,y)<=eps);
el2i=[el2 eli(1,:)];
j=floor((el2i-1)/n);
i=el2i-1-j*n;
x=[];
y=[];
for l=1:p+1
    for m=1:p+1
        x=[x;(xn(m)+i)*h];
        y=[y;(xn(l)+j)*h];
    end
end
uh2=u([ltgmap2 ltgmapi2]);
u2=ue2(x,y);
ui([ltgmap2 ltgmapi2])=u2;
uh2=uh2.*(lvsfun(x,y)>=-eps);%Get rid of DOFs in Omega_2.
u2=u2.*(lvsfun(x,y)>=-eps);
%err=max(max(max(abs(uh1-u1)))/max(max(abs(u1))),max(max(abs(uh2-u2)))/max(max(abs(u2))));
errinf=max(max(abs([uh1-u1 uh2-u2])))/max(max(abs([u1 u2])));




% Relative error in L^2 norm
err=0;
erri=0; % Interpolation error
uL2=0;
for j=1:2
    ltgmap=zeros((p+1)^2,n^2);
    if j==1
        el=el1;
        ltgmap(:,el)=ltgmap1;
        ltgmap(:,eli(1,:))=ltgmapi1;
        ue=ue1;
        sgn=-1;
    else
        el=el2;
        ltgmap(:,el)=ltgmap2;
        ltgmap(:,eli(1,:))=ltgmapi2;
        ue=ue2;
        sgn=1;
    end
    jl=floor((el-1)/n);
    il=el-1-jl*n;
    x0=il.'*h; 
    x1=x0+h;
    y0=jl.'*h;
    y1=y0+h;
    fun=@(x,y) abs(ue(x,y)-upfemeleevaluh(x,y,u,ltgmap,n,p,xn)).^2;
    funi=@(x,y) abs(ue(x,y)-upfemeleevaluh(x,y,ui,ltgmap,n,p,xn)).^2;
    q=quad2(fun,x0,x1,y0,y1,p);
    err=err+sum(q);
    q=quad2(funi,x0,x1,y0,y1,p);
    erri=erri+sum(q);
    q=quad2(@(x,y) abs(ue(x,y)).^2,x0,x1,y0,y1,p);
    uL2=uL2+sum(q);
    q=upfemitfquad2d(fun,lvsfun,lvsfunx,lvsfuny,eli,n,p,sgn);
    err=err+sum(q);
    q=upfemitfquad2d(funi,lvsfun,lvsfunx,lvsfuny,eli,n,p,sgn);
    erri=erri+sum(q);
    q=upfemitfquad2d(@(x,y) abs(ue(x,y)).^2,lvsfun,lvsfunx,lvsfuny,eli,n,p,sgn);
    uL2=uL2+sum(q);
end
err=sqrt(err);
erri=sqrt(erri);
uL2=sqrt(uL2);
err0=err/uL2;
erri0=erri/uL2;
 
 

ee=zeros(1,n^2);uu=zeros(1,n^2); 
    
    

% Relative error in energy norm
err=0;
erri=0; % Interpolation error
uHe=0;
for j=1:2
    ltgmap=zeros((p+1)^2,n^2);
    if j==1
        el=el1;
        ltgmap(:,el)=ltgmap1;
        ltgmap(:,eli(1,:))=ltgmapi1;
        uex=ue1x;
        uey=ue1y;
        sgn=-1;
    else
        el=el2;
        ltgmap(:,el)=ltgmap2;
        ltgmap(:,eli(1,:))=ltgmapi2;
        uex=ue2x;
        uey=ue2y;
        sgn=1;
    end
    jl=floor((el-1)/n);
    il=el-1-jl*n;
    x0=il.'*h; 
    x1=x0+h;
    y0=jl.'*h;
    y1=y0+h;
    funx=@(x,y) a(j)*abs(uex(x,y)-upfemeleevaluhx(x,y,u,ltgmap,n,p,xn)).^2;
    funy=@(x,y) a(j)*abs(uey(x,y)-upfemeleevaluhy(x,y,u,ltgmap,n,p,xn)).^2;
    funix=@(x,y) a(j)*abs(uex(x,y)-upfemeleevaluhx(x,y,ui,ltgmap,n,p,xn)).^2;
    funiy=@(x,y) a(j)*abs(uey(x,y)-upfemeleevaluhy(x,y,ui,ltgmap,n,p,xn)).^2;
    q=quad2(funx,x0,x1,y0,y1,p)+quad2(funy,x0,x1,y0,y1,p);
    err=err+sum(q);
    ee(el)=q;
    q=quad2(funix,x0,x1,y0,y1,p)+quad2(funiy,x0,x1,y0,y1,p);
    erri=erri+sum(q);
    q=quad2(@(x,y) a(j)*(abs(uex(x,y)).^2+abs(uey(x,y)).^2),x0,x1,y0,y1,p);
    uHe=uHe+sum(q);
    uu(el)=q;
    q= upfemitfquad2d(funx,lvsfun,lvsfunx,lvsfuny,eli,n,p,sgn)...
      +upfemitfquad2d(funy,lvsfun,lvsfunx,lvsfuny,eli,n,p,sgn);
    err=err+sum(q);
    ee(eli(1,:))=ee(eli(1,:))+q;
    q= upfemitfquad2d(funix,lvsfun,lvsfunx,lvsfuny,eli,n,p,sgn)...
      +upfemitfquad2d(funiy,lvsfun,lvsfunx,lvsfuny,eli,n,p,sgn);
    erri=erri+sum(q);
    q=upfemitfquad2d(@(x,y) a(j)*(abs(uex(x,y)).^2+abs(uey(x,y)).^2),lvsfun,lvsfunx,lvsfuny,eli,n,p,sgn);
    uHe=uHe+sum(q);
    uu(eli(1,:))=uu(eli(1,:))+q;
end
err=sqrt(err);
erri=sqrt(erri);
uHe=sqrt(uHe);
err1=err/uHe;
erri1=erri/uHe;


%  figure;
%  set(gcf,'position',[200,200,400,400])
%     ee=sqrt(ee);
%     uu=sqrt(uu);
%     ee=ee./uu;
%     ee=reshape(ee,n,n);
%     ee_pad = nan(n+1, n+1); 
%     ee_pad(1:end-1, 1:end-1) = ee; 
%     eex=linspace(0,1,n+1);
%     [xx1,yy1]=meshgrid(eex,eex);
%     pcolor(xx1,yy1,ee_pad');
%     cb = colorbar;       % 画出色标
% cb.FontSize = 15;
%     axis off;
%     axis tight;
%     axis equal;
%     hold on;
%     th=0:0.01:2*pi;
% % x=1/2+cos(th)/sqrt(9);
% % y=1/2+sin(th)/sqrt(9);
% x=1/2+(1/2+sin(5*th)/7).*cos(th)/2;
% y=1/2+(1/2+sin(5*th)/7).*sin(th)/2;
% plot(x,y,'r','LineWidth',1);
% hold on;
% k_val = 0:4;
% th_max_kappa = (3*pi/10) + k_val*2*pi/5;
% r_min = 1/2 * (1/2 - 1/7);
% x_max = 1/2 + r_min * cos(th_max_kappa);
% y_max = 1/2 + r_min * sin(th_max_kappa);
% plot(x_max, y_max, 'co', 'MarkerSize', 4, 'MarkerFaceColor','c'); % 红色圆圈标出极大曲率点
% hold off;

% if normflag==2
%     errx=0;
%     erry=0;
%     uHe=0;
%     for m=1:length(el1)
%         em=el1(m);
%         jm=floor((em-1)/n);
%         im=em-1-jm*n;
%         funx=@(x,y) a(1)^2*abs(ue1x(x,y)-upfemeleevaluhx(x,y,u,em,ltgmap1(:,m),n,p,xn)).^2;
%         funy=@(x,y) a(1)^2*abs(ue1y(x,y)-upfemeleevaluhy(x,y,u,em,ltgmap1(:,m),n,p,xn)).^2;
%         errx=errx+quad2(@(x,y) funx((im+x)*h,(jm+y)*h),0,1,0,1,p);
%         erry=erry+quad2(@(x,y) funy((im+x)*h,(jm+y)*h),0,1,0,1,p); 
%         uHe=uHe+quad2(@(x,y) a(1)^2*abs(ue1x((im+x)*h,(jm+y)*h)).^2,0,1,0,1,p);
%         uHe=uHe+quad2(@(x,y) a(1)^2*abs(ue1y((im+x)*h,(jm+y)*h)).^2,0,1,0,1,p);
%     end
%     for m=1:length(el2)
%         em=el2(m);
%         jm=floor((em-1)/n);
%         im=em-1-jm*n;
%         funx=@(x,y) a(2)^2*abs(ue2x(x,y)-upfemeleevaluhx(x,y,u,em,ltgmap2(:,m),n,p,xn)).^2;
%         funy=@(x,y) a(2)^2*abs(ue2y(x,y)-upfemeleevaluhy(x,y,u,em,ltgmap2(:,m),n,p,xn)).^2;
%         errx=errx+quad2(@(x,y) funx((im+x)*h,(jm+y)*h),0,1,0,1,p);  
%         erry=erry+quad2(@(x,y) funy((im+x)*h,(jm+y)*h),0,1,0,1,p); 
%         uHe=uHe+quad2(@(x,y) a(2)^2*abs(ue2x((im+x)*h,(jm+y)*h)).^2,0,1,0,1,p); 
%         uHe=uHe+quad2(@(x,y) a(2)^2*abs(ue2y((im+x)*h,(jm+y)*h)).^2,0,1,0,1,p);
%     end
%     for m=1:length(eli)
%         em=eli(1,m);
%         elitype=eli(2,m);
%         jm=floor((em-1)/n);
%         im=em-1-jm*n;
%         reflvsfun=@(x,y) lvsfun((im+x)*h,(jm+y)*h);
%         reflvsfunx=@(x,y) h*lvsfunx((im+x)*h,(jm+y)*h); 
%         reflvsfuny=@(x,y) h*lvsfuny((im+x)*h,(jm+y)*h);
%         funx=@(x,y) a(1)^2*abs(ue1x(x,y)-upfemeleevaluhx(x,y,u,em,ltgmapi1(:,m),n,p,xn)).^2;
%         funy=@(x,y) a(1)^2*abs(ue1y(x,y)-upfemeleevaluhy(x,y,u,em,ltgmapi1(:,m),n,p,xn)).^2;
%         errx=errx+upfemitfquad2d(@(x,y) funx((im+x)*h,(jm+y)*h),reflvsfun,reflvsfunx,reflvsfuny,elitype,p,-1);
%         erry=erry+upfemitfquad2d(@(x,y) funy((im+x)*h,(jm+y)*h),reflvsfun,reflvsfunx,reflvsfuny,elitype,p,-1);
%         funx=@(x,y) a(2)^2*abs(ue2x(x,y)-upfemeleevaluhx(x,y,u,em,ltgmapi2(:,m),n,p,xn)).^2;
%         funy=@(x,y) a(2)^2*abs(ue2y(x,y)-upfemeleevaluhy(x,y,u,em,ltgmapi2(:,m),n,p,xn)).^2;
%         errx=errx+upfemitfquad2d(@(x,y) funx((im+x)*h,(jm+y)*h),reflvsfun,reflvsfunx,reflvsfuny,elitype,p,1);
%         erry=erry+upfemitfquad2d(@(x,y) funy((im+x)*h,(jm+y)*h),reflvsfun,reflvsfunx,reflvsfuny,elitype,p,1);
%         uHe=uHe+upfemitfquad2d(@(x,y) a(1)^2*abs(ue1x((im+x)*h,(jm+y)*h)).^2,reflvsfun,reflvsfunx,reflvsfuny,elitype,p,-1);
%         uHe=uHe+upfemitfquad2d(@(x,y) a(1)^2*abs(ue1y((im+x)*h,(jm+y)*h)).^2,reflvsfun,reflvsfunx,reflvsfuny,elitype,p,-1);
%         uHe=uHe+upfemitfquad2d(@(x,y) a(2)^2*abs(ue2x((im+x)*h,(jm+y)*h)).^2,reflvsfun,reflvsfunx,reflvsfuny,elitype,p,1);
%         uHe=uHe+upfemitfquad2d(@(x,y) a(2)^2*abs(ue2y((im+x)*h,(jm+y)*h)).^2,reflvsfun,reflvsfunx,reflvsfuny,elitype,p,1);
%     end
%     
%     uHe=sqrt(uHe);
%     err=sqrt(errx+erry);
%     err=err/uHe;
% 
% end


% Relative error in flux
err=0;
erri=0; % Interpolation error
uHe=0;
f0=0;
for j=1:2
    ltgmap=zeros((p+1)^2,n^2);
    if j==1
        el=el1;
        ltgmap(:,el)=ltgmap1;
        ltgmap(:,eli(1,:))=ltgmapi1;
        uex=ue1x;
        uey=ue1y;
        sgn=-1;
        fj=f1;
    else
        el=el2;
        ltgmap(:,el)=ltgmap2;
        ltgmap(:,eli(1,:))=ltgmapi2;
        uex=ue2x;
        uey=ue2y;
        sgn=1;
        fj=f2;
    end
    jl=floor((el-1)/n);
    il=el-1-jl*n;
    x0=il.'*h; 
    x1=x0+h;
    y0=jl.'*h;
    y1=y0+h;
    funx=@(x,y) a(j)^2*abs(uex(x,y)-upfemeleevaluhx(x,y,u,ltgmap,n,p,xn)).^2;
    funy=@(x,y) a(j)^2*abs(uey(x,y)-upfemeleevaluhy(x,y,u,ltgmap,n,p,xn)).^2;
    funix=@(x,y) a(j)^2*abs(uex(x,y)-upfemeleevaluhx(x,y,ui,ltgmap,n,p,xn)).^2;
    funiy=@(x,y) a(j)^2*abs(uey(x,y)-upfemeleevaluhy(x,y,ui,ltgmap,n,p,xn)).^2;
    q=quad2(funx,x0,x1,y0,y1,p)+quad2(funy,x0,x1,y0,y1,p);
    err=err+sum(q);
    q=quad2(funix,x0,x1,y0,y1,p)+quad2(funiy,x0,x1,y0,y1,p);
    erri=erri+sum(q);
    q=quad2(@(x,y) a(j)^2*(abs(uex(x,y)).^2+abs(uey(x,y)).^2),x0,x1,y0,y1,p);
    uHe=uHe+sum(q);
    q=quad2(@(x,y) fj(x,y).^2,x0,x1,y0,y1,p);
    f0=f0+sum(q);
    q= upfemitfquad2d(funx,lvsfun,lvsfunx,lvsfuny,eli,n,p,sgn)...
      +upfemitfquad2d(funy,lvsfun,lvsfunx,lvsfuny,eli,n,p,sgn);
    err=err+sum(q);
    q= upfemitfquad2d(funix,lvsfun,lvsfunx,lvsfuny,eli,n,p,sgn)...
      +upfemitfquad2d(funiy,lvsfun,lvsfunx,lvsfuny,eli,n,p,sgn);
    erri=erri+sum(q);
    q=upfemitfquad2d(@(x,y) a(j)^2*(abs(uex(x,y)).^2+abs(uey(x,y)).^2),lvsfun,lvsfunx,lvsfuny,eli,n,p,sgn);
    uHe=uHe+sum(q);
    q=upfemitfquad2d(@(x,y) fj(x,y).^2,lvsfun,lvsfunx,lvsfuny,eli,n,p,sgn);
    f0=f0+sum(q);
end
err=sqrt(err);
erri=sqrt(erri);
uHe=sqrt(uHe);
% f0=sqrt(f0);

f0=3.434905181322524;
f1=24.494873722309165;
f2=1.479693423860519e+02;
if p==1
    errflux=err/f0;
elseif p==2
    errflux=err/f1;
elseif p==3
    errflux=err/f2;
end
%errflux=err/(f0+gbd0);
%errif=erri/uHe;


