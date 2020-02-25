function [U,Q,detA] = Hholder(A)
%Householder Calcul d'une décomposition QU Householder de A
%   Q orthogonale, U triangulaire supérieure
[n,m] = size(A);
if(n ~= m)
    error("Valable uniquement pour les matrices carrées");
end
clear m;

Q = eye(n);
signat = 0; % Pour le calcul du déterminant de 'Q'

U = A;

for i = 1:(n-1)
    
   % S'arrête si la matrice U est triangulaire supérieure
   if istriu(U)
       break;
   end
   
   % v vecteur de la matrice de Householder
   v = zeros(n,1);
   
   p = norm(U(i:n,i));
   v(i) = U(i,i)+sign(U(i,i))*p;
   v((i+1):n) = U((i+1):n,i);
   
   
   % Calcul de la norme inclus dans H pour éviter les erreurs d'arrondi
   H = eye(n) - (2*(v*v')) / (norm(v)^2);
   signat = signat + 1;
   
   Q = Q*H;
   U = H*U;
end

detA = (-1).^(mod(signat,2)) * prod(diag(U));
if(abs(detA) > eps)
       disp("La matrice A est inversible : le produit est donc unique.");
       disp("(En considérant une diagonale de U positive)");

% Pour terminer on modifie U de manière à avoir uniquement des nombres
% positifs sur la diagonale.
   P = diag(diag(U) ./ abs(diag(U)));
   U = P*U;
   Q = Q*P;
    
else
    disp("La matrice n'est pas inversible : la décomposition n'est donc pas forcément unique !");
    detA = 0;
end
end