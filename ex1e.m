%Infulunce of jump of coefficient
A=[1e6 1;1e5 1;1e4 1;1e3 1;1e2 1;10 1;1 1;1 10;1 1e2;1 1e3;1 1e4;1 1e5;1 1e6;1 1e7;1 1e8];
na=size(A,1);
theta=[0 0];
P=1:3;
n=32;
beta=1;
gamma0=100;
gamma1=100;
fname=['ex1data_' num2str(beta) '_' num2str(gamma0) '_n' num2str(n)];
for p=P
    fprintf('\n p=%g,  n=%g.\n',p,n);
    Err1=zeros(1,na);
    Err0=zeros(1,na);
    Errinf=zeros(1,na);
    Errflux=zeros(1,na);
    Cond=zeros(2,na);
    for j=1:na
        a=A(j,:);
        fprintf('  a1=%g, a2=%g\n',a(1),a(2));
        %[err1,err0,errinf,u,el1,el2,eli,ltg,cn,en1_2,en1_4,en2_2,en2_4]=ex1fun(a,n,p,beta,gamma0,gamma1);
        [err1,err0,errinf,errflux,u,el1,el2,eli,ltg,cn,en1_2,en1_4,en2_2,en2_4]=ex1fun(a,theta,n,p,beta,gamma0,gamma1);
        Err1(j)=err1;
        Err0(j)=err0;
        Errinf(j)=errinf;
        Errflux(j)=errflux;
        Cond(:,j)=cn;
    end
    save([fname '_p' num2str(p) '.mat'],'A','n','p','beta','gamma0','Err1','Err0','Errinf','Cond','Errflux');
end