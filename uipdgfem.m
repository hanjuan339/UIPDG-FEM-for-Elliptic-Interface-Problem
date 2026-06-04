function [u,ell1,ell2,eli,ltg,cn,en1_2,en1_4,en2_2,en2_4,e,ltgmaps]=uipdgfem(a,theta,f1,f2,gd,gn,lvsfun,lvsfunx,lvsfuny,n,p,beta,gamma0,gamma1)
% domain:[0,1]*[0,1];


h=1/n;

% 排列整体节点编号。
[ell1,ell2,eli]=upfemeleide(lvsfun,lvsfunx,lvsfuny,n);
[en1_2,en1_4,en2_2,en2_4,e]=upfemcb_ele(eli,lvsfun,lvsfunx,lvsfuny,n,p);
el1=1:n^2;
el2=1:n^2;
el1([ell2 en1_2(:)' en1_4(:)'])=[];
el2([ell1 en2_2(:)' en2_4(:)'])=[];

[ltgmaps,ltg,ndofs,sn1,sn2]=upfemdofs(ell1,ell2,eli,el1,el2,en1_2,en1_4,en2_2,en2_4,n,p);


% 标记需要加罚的内部边。
%e=mark_e(en1,en2,lvsfun,n);
% [E1,E2,Ei1,Ei2]=mark_e(n,en1_2,en1_4,en2_2,en2_4,el1,el2,ell1,ell2,eli);

[K1,K2,Ki11,Ki12,Ki21,Ki22,MM,M12,M11,J12,J21,J11,J22]=upfemasmcK(a,theta,lvsfun,lvsfunx,lvsfuny,ell1,ell2,eli,e,n,p,beta,gamma0,gamma1);

ldofs=(p+1)^2;

K=sparse(ndofs,ndofs);
for l=1:ldofs
    for k=1:ldofs
        m=k+(l-1)*ldofs;
        K=K+sparse(ltg{1}(l,:),ltg{1}(k,:),K1(m,:),ndofs,ndofs);
        K=K+sparse(ltg{2}(l,:),ltg{2}(k,:),K2(m,:),ndofs,ndofs);
        K=K+sparse(ltg{3}(l,:),ltg{3}(k,:),Ki11(m,:),ndofs,ndofs);
        K=K+sparse(ltg{4}(l,:),ltg{4}(k,:),Ki22(m,:),ndofs,ndofs);
        K=K+sparse(ltg{4}(l,:),ltg{3}(k,:),Ki12(m,:),ndofs,ndofs);
        K=K+sparse(ltg{3}(l,:),ltg{4}(k,:),Ki21(m,:),ndofs,ndofs);
    end
end

in1=[-n 1 n -1];
in2=cell(1,4);
in2{1}=[1:p+1;p*(p+1)+1:(p+1)^2];
in2{2}=[p+1:p+1:(p+1)^2;1:p+1:p*(p+1)+1];
in2{3}=[in2{1}(2,:);in2{1}(1,:)];
in2{4}=[in2{2}(2,:);in2{2}(1,:)];
el={[ell1,eli(1,:)],[ell2,eli(1,:)]};
Ltg={[ltg{1} ltg{3}],[ltg{2} ltg{4}]};
in3=[3 4 1 2];

Ii=repmat((1:ldofs)',1,ldofs);
Ii=Ii(:);
Ji=repmat(1:ldofs,ldofs,1);
Ji=Ji(:);


for l=1:2
    for i=1:size(e{2*l-1},2)
        for j=1:4
            if e{2*l-1}(j+1,i)~=0
                ii=ltgmaps{2+l}(in2{j}(1,:),i);
                ij=Ltg{l}(in2{j}(2,:),el{l}==e{2*l-1}(1,i)+in1(j));
                II=repmat(ii,1,p+1);
                IJ=repmat(ij',p+1,1);
                K=K+sparse(IJ(:),II(:),-MM(:)*a(l),ndofs,ndofs);
                II=repmat(ij,1,p+1);
                IJ=repmat(ii',p+1,1);
                K=K+sparse(IJ(:),II(:),-MM(:)*a(l),ndofs,ndofs);
                II=repmat(ii,1,p+1);
                IJ=repmat(ii',p+1,1);
                K=K+sparse(IJ(:),II(:),MM(:)*a(l),ndofs,ndofs);
                II=repmat(ij,1,p+1);
                IJ=repmat(ij',p+1,1);
                K=K+sparse(IJ(:),II(:),MM(:)*a(l),ndofs,ndofs);
                ei=el{l}==e{2*l-1}(1,i)+in1(j);
                K=K+sparse(Ltg{l}(Ji,ei),ltgmaps{2+l}(Ii,i),-a(l)*M12{j}(:),ndofs,ndofs);
                K=K+sparse(ltgmaps{2+l}(Ji,i),Ltg{l}(Ii,ei),-a(l)*M12{in3(j)}(:),ndofs,ndofs);
                K=K+sparse(ltgmaps{2+l}(Ji,i),ltgmaps{2+l}(Ii,i),-a(l)*M11{j}(:),ndofs,ndofs);
                K=K+sparse(Ltg{l}(Ji,ei),Ltg{l}(Ii,ei),-a(l)*M11{in3(j)}(:),ndofs,ndofs);
            end
             if e{2*l}(j+1,i)~=0
                ej=eli(1,:)==e{2*l}(1,i)+in1(j);
                K=K+sparse(ltg{2+l}(Ji,ej),ltgmaps{2+l}(Ii,i),J12{l}(:,4*(i-1)+j),ndofs,ndofs);
                K=K+sparse(ltgmaps{2+l}(Ji,i),ltg{2+l}(Ii,ej),J21{l}(:,4*(i-1)+j),ndofs,ndofs);
                K=K+sparse(ltgmaps{2+l}(Ji,i),ltgmaps{2+l}(Ii,i),J11{l}(:,4*(i-1)+j),ndofs,ndofs);
                K=K+sparse(ltg{2+l}(Ji,ej),ltg{2+l}(Ii,ej),J22{l}(:,4*(i-1)+j),ndofs,ndofs);
            end
        end
    end
end

if beta==1
    K=(K+K.')/2;
end


[F1,F2,Fi1,Fi2]=upfemasmfK(a,f1,f2,gd,gn,lvsfun,lvsfunx,lvsfuny,n,p,beta,gamma0);
F=sparse(ndofs,1);
for l=1:ldofs
    F=F+sparse(ltg{1}(l,:),1,F1(l,:),ndofs,1);
    F=F+sparse(ltg{2}(l,:),1,F2(l,:),ndofs,1);
    F=F+sparse(ltg{3}(l,:),1,Fi1(l,:),ndofs,1);
    F=F+sparse(ltg{4}(l,:),1,Fi2(l,:),ndofs,1);
end

%Homogeneous Dirichlet boundary condition
jj=floor((ell2-1)/n);
ii=ell2-1-jj*n;
x1=ii*h;
y1=jj*h;
x2=x1+h;
y2=y1+h;
nn=zeros(1,ndofs);
m= abs(x1)<=eps; %el2(m): left boundary elements
kk=ltg{2}(1:p+1:ldofs,m);
nn(kk)=1;
m= abs(x2-1)<=eps; %el2(m): right boundary elements
kk=ltg{2}(p+1:p+1:ldofs,m);
nn(kk)=1;
m= abs(y1)<=eps; %el2(m): lower boundary elements
kk=ltg{2}(1:p+1,m);
nn(kk)=1;
m= abs(y2-1)<=eps; %el2(m): upper boundary elements
kk=ltg{2}((p+1)*p+1:ldofs,m);
nn(kk)=1;

jj=floor((eli(1,:)-1)/n);
ii=eli(1,:)-1-jj*n;
x1=ii*h;
y1=jj*h;
x2=x1+h;
y2=y1+h;
m= abs(x1)<=eps; %el2(m): left boundary elements
kk=ltg{4}(1:p+1:ldofs,m);
nn(kk)=1;
m= abs(x2-1)<=eps; %el2(m): right boundary elements
kk=ltg{4}(p+1:p+1:ldofs,m);
nn(kk)=1;
m= abs(y1)<=eps; %el2(m): lower boundary elements
kk=ltg{4}(1:p+1,m);
nn(kk)=1;
m= abs(y2-1)<=eps; %el2(m): upper boundary elements
kk=ltg{4}((p+1)*p+1:ldofs,m);
nn(kk)=1;

bdofs= find(nn); %Boundary DOFs

B=speye(ndofs,ndofs);

%Remove slave nodes

for m=1:size(sn1,2)
    j=sn1(1,m);
    for k=2:p+2
        l=sn1(k,m);
        B(j,l)=B(j,l)+sn1(p+k+1,m)*B(j,j);
    end
end

for m=1:size(sn2,2)
    j=sn2(1,m);
    for k=2:(p+1)^2+1
        l=sn2(k,m);
        B(j,l)=B(j,l)+sn2((p+1)^2+k,m)*B(j,j);
    end
end

            
% B(:,bdofs)=[];

B(:,[sn1(1,:) sn2(1,:) bdofs])=[];

%Eliminate the Dirichlet DOFs
K=B.'*K*B;
cn(1)=condest(K);
F=B.'*F;
%Solve
fprintf('Solving the linear system...');
%Diagonal scaling 
d=diag(K);
D=sparse(1:length(d),1:length(d),sqrt(1./d));
K=D*K*D;
cn(2)=condest(K);
if beta==1
    K=(K+K.')/2;
end
F=D*F;
u=K\F;
u=D*u;
%Expand
u=B*u;
fprintf(' Done\n');



