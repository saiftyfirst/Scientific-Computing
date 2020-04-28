C = [2 5 3 8 7 6 1 4];
k = log2(length(C));
% for l = 1:k
%     for j = 1:(2^(k-l))
%         C(2^(l) * j) = C(2^(l) * j) + C(2^(l) * j + 2^(l-1));
%     end
% end
A = [1 2 3 4 5 6 7 8];
for i=2:8 
    A(i) = A(i-1) * A(i);
end
A


