function [en1_2,en1_4,en2_2,en2_4,e]=upfemcb_ele(eli,lvsfun,lvsfunx,lvsfuny,n,p)
%Combining some interface elements with theirs neighbours so that the area
%ratios are not small.

r0=1/4; %The minimum ratio tolerance r0 should >=1/3.
h=1/n; %Mesh size
% ldofs=(p+1)^2;
% Ltgmap1=zeros(ldofs,n^3);
% Ltgmap2=Ltgmap1;
% Ltgmap1(:,el1)=ltgmap1;
% Ltgmap2(:,el2)=ltgmap2;
% ei=eli(1,:);
% Ltgmap1(:,ei)=ltgmapi1;
% Ltgmap2(:,ei)=ltgmapi2;


%Computing aera ratios K_1/K
ri1=upfemitfquad2d(@(x,y) ones(size(x)),lvsfun,lvsfunx,lvsfuny,eli,n,p,-1)/h^2;

j=find(ri1<r0); %The parts in Om_1 are small
en1=zeros(2,length(j));
for m=1:length(j)
    em=eli(1,j(m));
    elitype=eli(2,j(m));
    jm=floor((em-1)/n);
    im=em-1-jm*n;
    x1=im*h;
    y1=jm*h;    
    %Find the neighbour element to combine
    switch elitype
        case 1
            emn=em+n;
        case 2
            emn=em-n;
        case 3
            emn=em+1;
        case 4
            emn=em-1;
        case 5
            if lvsfun(x1+h/2,y1+3*h/2)<lvsfun(x1+3*h/2,y1+h/2)
                emn=em+n;
            else
                emn=em+1;
            end
        case 6
            if lvsfun(x1+h/2,y1+3*h/2)<lvsfun(x1-h/2,y1+h/2)
                emn=em+n;
            else
                emn=em-1;
            end
        case 7
            if lvsfun(x1+h/2,y1-h/2)<lvsfun(x1+3*h/2,y1+h/2)
                emn=em-n;
            else
                emn=em+1;
            end
        case 8
            if lvsfun(x1+h/2,y1-h/2)<lvsfun(x1-h/2,y1+h/2)
                emn=em-n;
            else
                emn=em-1;
            end
        case {13,24}
            emn=em+n;
        case {15,22}
            emn=em-n;
        case {17,28}
            emn=em+1;
        case {19,26}
            emn=em-1;
    end
    en1(:,m)=[em;emn];
end

    

d1=unique(en1(:));
c1=histc(en1(:),d1);
b1=c1==2;
a1=d1(b1);
en1_4=zeros(1,length(a1));
for i=1:length(a1)
    I1=en1(1,:)==a1(i);
    I2=en1(2,:)==a1(i);
    if nnz(I1)
        en1_4(i)=min([a1(i) en1(1,I2) en1(2,I1) en1(1,I2)+en1(2,I1)-a1(i)]);
        en1(:,I1|I2|(en1(1,:)==en1(1,I2)+en1(2,I1)-a1(i)))=[];
    elseif nnz(I2)
        if sum(en1(1,I2))-sum(en1(2,I2))==0
            ele=sum(en1(2,I2))/2;
            elej=floor((ele-1)/n);
            elei=ele-1-elej*n;
            if sum(abs(en1(1,I2)-en1(2,I2)))==2
                if lvsfun(elei*h+h/2,elej*h+3*h/2)>lvsfun(elei*h+h/2,elej*h-h/2)
                    en1(2,I2)=en1(1,I2)-n;
                else
                    en1(2,I2)=en1(1,I2)+n;
                end
            else 
                if lvsfun(elei*h-h/2,elej*h+h/2)>lvsfun(elei*h+3*h/2,elej*h+h/2)
                    en1(2,I2)=en1(1,I2)+1;
                else
                    en1(2,I2)=en1(1,I2)-1;
                end
            end
            en1_4(i)=0;
        else
            en1_4(i)=min([a1(i) en1(1,I2) sum(en1(1,I2))-a1(i)]); 
            en1(:,I2|(en1(1,:)==sum(en1(1,I2))-a1(i)))=[];
        end
    else
        en1_4(i)=0;
    end
end
en1_2=en1;
en1_4=en1_4(en1_4~=0);
en1_4=[en1_4;en1_4+1;en1_4+n;en1_4+n+1];
            
j=find(ri1>1-r0); %The parts in Om_2 are small
en2=zeros(2,length(j));
for m=1:length(j)
    em=eli(1,j(m));
    elitype=eli(2,j(m));
    jm=floor((em-1)/n);
    im=em-1-jm*n;
    x1=im*h;
    y1=jm*h;    
    %Find the neighbour element to combine
    switch elitype
        case 1
            emn=em-n;
        case 2
            emn=em+n;
        case 3
            emn=em-1;
        case 4
            emn=em+1;
        case 9
            if lvsfun(x1+h/2,y1-h/2)>lvsfun(x1-h/2,y1+h/2)
                emn=em-n;
            else
                emn=em-1;
            end             
        case 10
            if lvsfun(x1+h/2,y1-h/2)>lvsfun(x1+3*h/2,y1+h/2)
                emn=em-n;
            else
                emn=em+1;
            end
        case 11
            if lvsfun(x1+h/2,y1+3*h/2)>lvsfun(x1-h/2,y1+h/2)
                emn=em+n;
            else
                emn=em-1;
            end
        case 12
            if lvsfun(x1+h/2,y1+3*h/2)>lvsfun(x1+3*h/2,y1+h/2)
                emn=em+n;
            else
                emn=em+1;
            end
        case {14,23}
            emn=em-n;
        case {16,21}
            emn=em+n;
        case {18,27}
            emn=em-1;
        case {20,25}
            emn=em+1;
    end
    en2(:,m)=[em;emn];
end

d1=unique(en2(:));
c1=histc(en2(:),d1);
b1=c1==2;
a1=d1(b1);
en2_4=zeros(1,length(a1));
for i=1:length(a1)
    I1=en2(1,:)==a1(i);
    I2=en2(2,:)==a1(i);
    if nnz(I1)
        en2_4(i)=min([a1(i) en2(1,I2) en2(2,I1) en2(1,I2)+en2(2,I1)-a1(i)]);
        en2(:,I1|I2|(en2(1,:)==en2(1,I2)+en2(2,I1)-a1(i)))=[];
    elseif nnz(I2)
        if sum(en2(1,I2))-sum(en2(2,I2))==0
            ele=sum(en2(2,I2))/2;
            elej=floor((ele-1)/n);
            elei=ele-1-elej*n;
            if sum(abs(en2(1,I2)-en2(2,I2)))==2
                if lvsfun(elei*h+h/2,elej*h+3*h/2)<lvsfun(elei*h+h/2,elej*h-h/2)
                    en2(2,I2)=en2(1,I2)-n;
                else
                    en2(2,I2)=en2(1,I2)+n;
                end
            else 
                if lvsfun(elei*h-h/2,elej*h+h/2)<lvsfun(elei*h+3*h/2,elej*h+h/2)
                    en2(2,I2)=en2(1,I2)+1;
                else
                    en2(2,I2)=en2(1,I2)-1;
                end
            end
            en2_4(i)=0;
        else
            en2_4(i)=min([a1(i) en2(1,I2) sum(en2(1,I2))-a1(i)]); 
            en2(:,I2|(en2(1,:)==sum(en2(1,I2))-a1(i)))=[];
        end
    else
        en2_4(i)=0;
    end
end
en2_2=en2;
en2_4=en2_4(en2_4~=0);
en2_4=[en2_4;en2_4+1;en2_4+n;en2_4+n+1];



%  ———
% |  3   |
% |4    2|  ：enl(i,k)单元四个边的排列.
% |  1   |
%  ———
e=cell(1,4);%四部分分别是区域1中完全没入边和部分没入边，以及区域2的相关边信息
            %e{1},e{3}是整边情况，e{2}，e{4}是相交边的情况
ell={en1_2,en2_2,en1_4,en2_4};
index=[-1 1];
for l=1:2
    e{2*l-1}=zeros(4,2*size(ell{l},2)+4*size(ell{2+l},2));
    e{2*l}=e{2*l-1};
    y1=floor(([ell{l}(:)' ell{l+2}(:)']-1)/n);
    x1=[ell{l}(:)' ell{l+2}(:)']-n*y1-1;
    x2=x1+1;y2=y1;
    x3=x2;y3=y2+1;
    x4=x1;y4=y1+1;
    d1=lvsfun(x1*h,y1*h)*index(l);
    d2=lvsfun(x2*h,y2*h)*index(l);
    d3=lvsfun(x3*h,y3*h)*index(l);
    d4=lvsfun(x4*h,y4*h)*index(l);
    k11=find(d1>=0&d2>=0);
    e{2*l-1}(1,k11)=ones(1,length(k11));
    k12=find((d1>0&d2<0)|(d1<0&d2>0));
    e{2*l}(1,k12)=ones(1,length(k12));
    k21=find(d2>=0&d3>=0);
    e{2*l-1}(2,k21)=ones(1,length(k21));
    k22=find((d2>0&d3<0)|(d2<0&d3>0));
    e{2*l}(2,k22)=ones(1,length(k22));
    k31=find(d3>=0&d4>=0);
    e{2*l-1}(3,k31)=ones(1,length(k31));
    k32=find((d3>0&d4<0)|(d3<0&d4>0));
    e{2*l}(3,k32)=ones(1,length(k32));
    k41=find(d4>=0&d1>=0);
    e{2*l-1}(4,k41)=ones(1,length(k41));
    k42=find((d4>0&d1<0)|(d4<0&d1>0));
    e{2*l}(4,k42)=ones(1,length(k42));
    
    k=ell{l}(1,:)-ell{l}(2,:);
    k1=find(k==1);
    e{2*l-1}(4,2*k1-1)=0;e{2*l}(4,2*k1-1)=0;
    e{2*l-1}(2,2*k1)=0;e{2*l}(2,2*k1)=0;
    k2=find(k==-1);
    e{2*l-1}(2,2*k2-1)=0;e{2*l}(2,2*k2-1)=0;
    e{2*l-1}(4,2*k2)=0;e{2*l}(4,2*k2)=0;
    k3=find(k==n);
    e{2*l-1}(1,2*k3-1)=0;e{2*l}(1,2*k3-1)=0;
    e{2*l-1}(3,2*k3)=0;e{2*l}(3,2*k3)=0;
    k4=find(k==-n);
    e{2*l-1}(3,2*k4-1)=0;e{2*l}(3,2*k4-1)=0;
    e{2*l-1}(1,2*k4)=0;e{2*l}(1,2*k4)=0;
    
    e{2*l-1}([2 3],(0:size(ell{l+2},2)-1)*4+2*size(ell{l},2)+1)=0;
    e{2*l}([2 3],(0:size(ell{l+2},2)-1)*4+2*size(ell{l},2)+1)=0;
    e{2*l-1}([3 4],(0:size(ell{l+2},2)-1)*4+2*size(ell{l},2)+2)=0;
    e{2*l}([3 4],(0:size(ell{l+2},2)-1)*4+2*size(ell{l},2)+2)=0;
    e{2*l-1}([1 2],(0:size(ell{l+2},2)-1)*4+2*size(ell{l},2)+3)=0;
    e{2*l}([1 2],(0:size(ell{l+2},2)-1)*4+2*size(ell{l},2)+3)=0;
    e{2*l-1}([1 4],(0:size(ell{l+2},2)-1)*4+2*size(ell{l},2)+4)=0;
    e{2*l}([1 4],(0:size(ell{l+2},2)-1)*4+2*size(ell{l},2)+4)=0;
    
    e{2*l-1}=[[ell{l}(:)' ell{l+2}(:)'];e{2*l-1}];
    e{2*l}=[[ell{l}(:)' ell{l+2}(:)'];e{2*l}];
end

in1=[-n 1 n -1];
in2=[3 4 1 2];
for l=1:2
    el=[e{2*l-1}(1,:);e{2*l-1}(2:5,:)+e{2*l}(2:5,:)];
    for i=2:size(el,2)
        for j=1:4
            if el(j+1,i)~=0
                a=find(el(1,1:i-1)==el(1,i)+in1(j),1);
                if ~isempty(a)
                    e{2*l-1}(1+in2(j),a)=0;
                    e{2*l}(1+in2(j),a)=0;
                end
            end
        end
    end
end

%去掉边界边
for l=1:2
    jj=floor((e{2*l-1}(1,:)-1)/n);
    ii=e{2*l-1}(1,:)-1-jj*n;
    x1=ii*h;
    y1=jj*h;
    x2=x1+h;
    y2=y1+h;
    m= find(abs(y1)<=eps); %el2(m): lower boundary elements
    e{2*l-1}(2,m)=zeros(size(m));
    m= find(abs(x2-1)<=eps); %el2(m): right boundary elements
    e{2*l-1}(3,m)=zeros(size(m));
    m= find(abs(y2-1)<=eps); %el2(m): upper boundary elements
    e{2*l-1}(4,m)=zeros(size(m));
    m= find(abs(x1)<=eps); %el2(m): left boundary elements
    e{2*l-1}(5,m)=zeros(size(m));
    
end