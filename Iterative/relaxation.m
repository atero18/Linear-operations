function [matIter,invM] = relaxation(A,omega)
[D,E,F] = partiesMat(A);
M = D ./ omega - E;
invM = inv(M);

matIter = invM * (((1 - omega) / omega) .* D +F);
end

