function [wOpti,rho] = wOptiHermPos(A)
% Renvoi la valeur optimale pour la relaxation de A
n = size(A);
n = n(1);
D = diag(partiesMat(A));
rho = max(abs(eig(eye(n) - diag(1./D)*A)));
wOpti = 2 / (1 + sqrt(1 - rho^2));
end

