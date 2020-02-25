function [matIter,invM] = Jacobi(A)
[D,E,F] = partiesMat(A);

invM = diag(1 ./ diag(D));
matIter = invM * (E + F);

end

