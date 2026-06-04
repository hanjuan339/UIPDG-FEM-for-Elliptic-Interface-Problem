function [ltgmaps,ltg,ndofs,sn1,sn2]=upfemdofs(ell1,ell2,eli,el1,el2,en1_2,en1_4,en2_2,en2_4,n,p)
%The map from the local DOFs to the global DOFs.
%   el1: elements inside the interface (where the lvsfun<0)
%   el2: elements outside the interface (where the lvsfun>0)
%   eli: interface elements. eli is a two-row matrix. The first row gives
%       the indices. The second row indicates which part the interface 
%       elements locate:
%   n: the unit square is divided into n x n subsquares of the same size.
%   p: the degree of element polynomials (in each x and y, respectively).
%   ltgmaps{1},ltgmaps{2}:el1 和 el2 中单元的局部编号和总体编号的一一对应；
%   ltgmaps{3},ltgmaps{4}:en1 和 en2 中单元的局部编号和总体编号的一一对应。
%   ndofs: number of DOFs

ldofs=(p+1)^2; %Number of DOFs in each element
ell={el1,el2,en1_2,en2_2,en1_4,en2_4};
ltgmaps=cell(1,6);
ndofs=0;
for l=1:2
    elli=ell{l}; %Elements in Omega_{lh}

    %Each DOF is located by (i,j): i-th column, j-th row.
    Il=zeros(ldofs,length(elli));
    Jl=zeros(ldofs,length(elli));
    jelli=floor((elli-1)/n);
    ielli=elli-1-jelli*n;
    for k=1:ldofs
        jk=floor((k-1)/(p+1))+1;
        ik=k-(jk-1)*(p+1);
        Il(k,:)=ielli*p+ik;
        Jl(k,:)=jelli*p+jk;
    end
    Cl=sparse(Il,Jl,1,n*ldofs,n*ldofs);
    [ii,jj]=find(Cl);
    Cl=sparse(ii,jj,ndofs+1:ndofs+length(ii),n*ldofs,n*ldofs);
    ndofs=ndofs+length(ii);
    ltgmaps{l}=zeros(size(Il));
    for k=1:ldofs
        ltgmaps{l}(k,:)=aij(Cl,Il(k,:),Jl(k,:));
    end
end


for l=3:6
    elli=ell{l};
    ltgmaps{l}=zeros(ldofs,size(elli,1)*size(elli,2));
    for i=1:size(elli,2)
        Il=zeros(ldofs,size(elli,1));
        Jl=zeros(ldofs,size(elli,1));
        jelli=floor((elli(:,i)-1)/n);
        ielli=elli(:,i)-1-jelli*n;
        for k=1:ldofs
            jk=floor((k-1)/(p+1))+1;
            ik=k-(jk-1)*(p+1);
            Il(k,:)=ielli*p+ik;
            Jl(k,:)=jelli*p+jk;
        end
        Cl=sparse(Il,Jl,1,n*ldofs,n*ldofs);
        [ii,jj]=find(Cl);
        Cl=sparse(ii,jj,ndofs+1:ndofs+length(ii),n*ldofs,n*ldofs);
        ndofs=ndofs+length(ii);
        for k=1:ldofs
            ltgmaps{l}(k,size(elli,1)*(i-1)+1:size(elli,1)*i)=aij(Cl,Il(k,:),Jl(k,:));
        end
    end
end

sn1=zeros(2*p+3,ndofs);
in_n=zeros(2*p+1,p+1);
inn=in_n;
in1=inn;
in_1=inn;
A=reshape(1:ldofs,p+1,p+1);
in_n(1:p+1,:)=A';
in_n(p+2:2*p+1,:)=A(:,2:p+1)';
inn(1:p+1,:)=flipud(A');
inn(p+2:2*p+1,:)=flipud(A(:,1:p)');
in1(1:p+1,:)=flipud(A);
in1(p+2:2*p+1,:)=flipud(A(1:p,:));
in_1(1:p+1,:)=A;
in_1(p+2:2*p+1,:)=A(2:p+1,:);

for l=1:2
    ee=ell{l+2};
    ei=ee(1,:)-ee(2,:);
    k=find(ei==n);
    aa1=zeros(2*p+1,(p+1)*length(k));
    aa1(1:p+1,:)=reshape(ltgmaps{l+2}(inn(1:p+1,:),2*k-1),p+1,(p+1)*length(k));
    aa1(p+2:2*p+1,:)=reshape(ltgmaps{l+2}(inn(p+2:2*p+1,:),2*k),p,(p+1)*length(k));
    k=find(ei==-n);
    aa2=zeros(2*p+1,(p+1)*length(k));
    aa2(1:p+1,:)=reshape(ltgmaps{l+2}(in_n(1:p+1,:),2*k-1),p+1,(p+1)*length(k));
    aa2(p+2:2*p+1,:)=reshape(ltgmaps{l+2}(in_n(p+2:2*p+1,:),2*k),p,(p+1)*length(k));
    k=find(ei==1);
    aa3=zeros(2*p+1,(p+1)*length(k));
    aa3(1:p+1,:)=reshape(ltgmaps{l+2}(in1(1:p+1,:),2*k-1),p+1,(p+1)*length(k));
    aa3(p+2:2*p+1,:)=reshape(ltgmaps{l+2}(in1(p+2:2*p+1,:),2*k),p,(p+1)*length(k));
    k=find(ei==-1);
    aa4=zeros(2*p+1,(p+1)*length(k));
    aa4(1:p+1,:)=reshape(ltgmaps{l+2}(in_1(1:p+1,:),2*k-1),p+1,(p+1)*length(k));
    aa4(p+2:2*p+1,:)=reshape(ltgmaps{l+2}(in_1(p+2:2*p+1,:),2*k),p,(p+1)*length(k));
    aa=[aa1 aa2 aa3 aa4];
    for j=1:p
        sn1(:,aa(2*j,:))=[aa(2*j,:);aa(1:2:2*p+1,:);repmat(coeff((2*j-1)/(2*p),p),1,size(aa,2))];
    end
end
sn1(:,sn1(1,:)==0)=[];

sn2=zeros(2*ldofs+1,ndofs);
aa=[ltgmaps{5}(1,1:4:4*size(en1_4,2)) ltgmaps{6}(1,1:4:4*size(en2_4,2))];
aa=repmat(aa,(2*p+1)^2,1)+repmat((0:(2*p+1)^2-1)',1,length(aa));
bb=(1:2:2*p+1)'+(0:2:2*p)*(2*p+1);
bb=bb(:);
cc=1:(2*p+1)^2;
cc(bb)=[];
ccy=floor((cc-1)/(2*p+1));
ccx=cc-ccy*(2*p+1)-1;
ccy=ccy/(2*p);
ccx=ccx/(2*p);
dd=zeros(ldofs,length(cc));
xn=linspace(0,1,p+1)';
for j=1:p+1
    for i=1:p+1
        dd(i+(j-1)*(p+1),:)=upfemrefbas(ccx,ccy,i,j,p,xn);
    end
end
for j=1:length(cc)
    sn2(:,aa(cc(j),:))=[aa(cc(j),:);aa(bb,:);repmat(dd(:,j),1,size(aa,2))];
end
sn2(:,sn2(1,:)==0)=[];


ltgmaps{3}=[ltgmaps{3} ltgmaps{5}];
ltgmaps{4}=[ltgmaps{4} ltgmaps{6}];
ltgmaps=ltgmaps(1:4);
aa=[en1_2(:)' en1_4(:)'];
bb=[en2_2(:)' en2_4(:)'];

ltg=cell(1,4);
ltg{1}=zeros(ldofs,length(ell1));% \Omega_1 内部
ltg{2}=zeros(ldofs,length(ell2));% \Omega_2 内部
ltg{3}=zeros(ldofs,size(eli,2));% \Omega_1 部分
ltg{4}=ltg{3};                 % \Omega_2 部分

[A,B]=find_matrix(el1,ell1);
[C,D]=find_matrix(ell1,el1(A));
E=find_matrix(aa,ell1(D));
[F,G]=find_matrix(eli(1,:),el1(B));
H=find_matrix(aa,eli(1,G));

ltg{1}(:,C)=ltgmaps{1}(:,A);
ltg{1}(:,D)=ltgmaps{3}(:,E);
ltg{3}(:,F)=ltgmaps{1}(:,B);
ltg{3}(:,G)=ltgmaps{3}(:,H);


[A,B]=find_matrix(el2,ell2);
[C,D]=find_matrix(ell2,el2(A));
E=find_matrix(bb,ell2(D));
[F,G]=find_matrix(eli(1,:),el2(B));
H=find_matrix(bb,eli(1,G));

ltg{2}(:,C)=ltgmaps{2}(:,A);
ltg{2}(:,D)=ltgmaps{4}(:,E);
ltg{4}(:,F)=ltgmaps{2}(:,B);
ltg{4}(:,G)=ltgmaps{4}(:,H);

