function cv=upfemcv(fun,A,B)
%Estimate the curvature of a function between A and B
%   A and B may be column vectors of the same size

if isempty(A)
    cv=[];
    return;
end

x1=A;
x2=B;
x3=(A+B)/2;
y1=fun(x1);
y2=fun(x2);
y3=fun(x3);
cv=2*abs((x2-x1).*(y3-y1)-(x3-x1).*(y2-y1))./sqrt(((x1-x2).^2+(y1-y2).^2).*((x2-x3).^2+(y2-y3).^2).*((x3-x1).^2+(y3-y1).^2));

% cv=2*abs((x2-x1).*(y3-y1)-(x3-x1).*(y2-y1)).*sqrt(((x1-x2).^2+(y1-y2).^2))./abs(x1-x2)...
%     ./sqrt(((x1-x2).^2+(y1-y2).^2).*((x2-x3).^2+(y2-y3).^2).*((x3-x1).^2+(y3-y1).^2));

