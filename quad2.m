function q = quad2(fun,A,B,c,d,p)
%Numerically evaluate double integral for p-th order UIPFEM.

%   Q = QUAD2(FUN,A,B,C,D,P) approximates the integral of
%   FUN(X,Y) over planar regions A <=X<= B and C(X) <=Y<=D(X), 
%   where A and B may be column vectors of the same size. 
%   FUN is a function handle, C and D may each be a column vector
%   or a function handle.

if isempty(A)
    q=[];
    return;
end

n=size(A,1);
K=2*ones(n,1);
if isa(c,'function_handle')
    cv=upfemcv(@(x) c(x,1:n),A,B); % Curvature between A and B
    K=max(K,ceil(cv/5)); %K+(cv>=30);
end
if isa(d,'function_handle')
    cv=upfemcv(@(x) d(x,1:n),A,B); % Curvature between A and B
    K=max(K,ceil(cv/5)); %K+(cv>=30);
end
q=zeros(n,1);
for k=1:max(K)
    l=(K==k);
    if nnz(l)
        m=k*p+1;%ceil(3*p/2);%k*p+1;
        [x,w]=lgwt(m,0,1);%Gaussian nodes and weights in [0,1]
        xn=repmat(A(l),1,m)+(B(l)-A(l))*x.';
        if isa(c,'function_handle')
            Cn=zeros(size(xn));
            for j=1:m
                Cn(:,j)=c(xn(:,j),l);
            end
        elseif isfloat(c) && size(c,1)==n
            Cn=repmat(c(l),1,m);
        else
            error(message('quad2:invalidC'));
        end
        if isa(d,'function_handle')
            Dn=zeros(size(xn));
            for j=1:m
                Dn(:,j)=d(xn(:,j),l);
            end
        elseif isfloat(d) && size(d,1)==n
            Dn=repmat(d(l),1,m);
        else
            error(message('quad2:invalidD'));
        end
        for j=1:m
            Y=repmat(Cn(:,j),1,m)+(Dn(:,j)-Cn(:,j))*x';%plot(repmat(xn(:,j),1,m),Y,'.');
            F=fun(repmat(xn(:,j),1,m),Y); %plot(Y,repmat(xn(:,j),1,m),'.');
            q(l)=q(l)+(w(j).*(Dn(:,j)-Cn(:,j))).*F*w;
        end
        q(l)=(B(l)-A(l)).*q(l);
    end
end
