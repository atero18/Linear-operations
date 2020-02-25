function [U,Q,detA] = Hholder(A)
%Householder Calcul d'une d�composition QU Householder de A
%   Q orthogonale, U triangulaire sup�rieure
[n,m] = size(A);
if(n ~= m)
    error("Valable uniquement pour les matrices carr�es");
end
clear m;

Q = eye(n);
signat = 0; % Pour le calcul du d�terminant de 'Q'

U = A;

for i = 1:(n-1)
    
   % S'arr�te si la matrice U est triangulaire sup�rieure
   if istriu(U)
       break;
   end
   
   % v vecteur de la matrice de Householder
   v = zeros(n,1);
   
   p = norm(U(i:n,i));
   v(i) = U(i,i)+sign(U(i,i))*p;
   v((i+1):n) = U((i+1):n,i);
   
   
   % Calcul de la norme inclus dans H pour �viter les erreurs d'arrondi
   H = eye(n) - (2*(v*v')) / (norm(v)^2);
   signat = signat + 1;
   
   Q = Q*H;
   U = H*U;
end

detA = (-1).^(mod(signat,2)) * prod(diag(U));
if(abs(detA) > eps)
       disp("La matrice A est inversible : le produit est donc unique.");
       disp("(En consid�rant une diagonale de U positive)");

% Pour terminer on modifie U de mani�re � avoir uniquement des nombres
% positifs sur la diagonale.
   P = diag(diag(U) ./ abs(diag(U)));
   U = P*U;
   Q = Q*P;
    
else
    disp("La matrice n'est pas inversible : la d�composition n'est donc pas forc�ment unique !");
    detA = 0;
end
end