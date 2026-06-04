function [A1,A2]=find_matrix(b,b1)
% b1 is a subset of b
% A is the index of b1 in b

A1=[];

for i=1:length(b1)
    if find(b==b1(i))
        A1=[A1 find(b==b1(i))];
    end
end

A2=1:length(b);
A2(A1)=[];
