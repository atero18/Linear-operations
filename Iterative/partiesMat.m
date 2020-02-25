function [D,E,F] = partiesMat(A)
% Renvoie la d�composition D,E,F de la matrice A
D = diag(diag(A));
E = -tril(A,-1);
F = -triu(A,1);
end

