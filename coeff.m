function coeff=coeff(x,p)
xn=linspace(0,1,p+1);
coeff=zeros(p+1,1);
for i=1:p+1
    xx=xn;
    xxi=repmat(x,1,p+1)-xn;
    xxi(i)=[];
    xx(i)=[];
    xxn=repmat(xn(i),1,p)-xx;
    coeff(i)=prod(xxi)./prod(xxn);
end
