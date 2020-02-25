function [matIter,invM] = Gauss_Seidel(A)
[D,E,F] = partiesMat(A);

M = D - E;
invM = inv(M);
matIter = invM * F;
end

