function [A,B] = creatAB(n)
    A = 2*diag(ones(n,1)) - diag(ones(n-1,1),-1) - diag(ones(n-1,1),1);
    B = zeros(1,n)';
    B([1,n]) = 1;
end

