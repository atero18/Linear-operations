function [U, L, P,Q, det] = gauss(A,pivotTech)
%pivot de Gauss : Rend une d�composition P'LUQ' de la matrice A ainsi que
%son d�terminant.
%   Par d�faut la m�thode est sur pivot partiel 'n�cessaire'.
% Pivot partiel et pivot total peuvent �tre forc�s.

[n,m] = size(A);
if n ~= m
    error("Valable uniquement pour les matrices carr�es");
end
clear m;

% M�thode pivot : 0 = pivot partiel (d�faut) ; 1 = pivot partiel forc� ;
% 2 = pivot total
if nargin == 1
    pivotTech = 0;
end

% Sert � conserver les changements de lignes / colonnes 
Q = eye(n);
signatureQ = 0;
P = eye(n);
signatureP = 0;

L = eye(n);
U = A;
needTotal = false;

for i = 1:(n-1)
    % Si le pivot est nul ou si le pivotage est forc� on le modifie
    if(pivotTech ~= 0 || abs(U(i,i)) < eps)
        l = i; %indice ligne pivot
        m = i; %indice colonne pivot (total)
        if(pivotTech <= 1)
            [M,l] = max(abs(U(i:n,i)));
            if(M == 0)
                needTotal = true;
            end
            % On change la position de la ligne pivot. Formule due � la
            % forme de la recherche (on prend le premier indice disponible)
            l = l(1) + i - 1;
        end
        if(needTotal || pivotTech == 2)
            needTotal = false;
            Val = abs(U(i:n,i:n));
            M = max(Val(:));
            
            if(M == 0)
                error("Aucun pivot disponible.");
            end
            [l,m] = find(abs(U(i:n,i:n)) == M);
            % Formule du � la forme de la liste du 'find'
            l = l(1) + i - 1;
            m = m(1) + i - 1;
            
        end
        
        % Cas o� ou intervertit 2 lignes
        if(l ~= i)
            signatureP = signatureP + 1;
            Plin = eye(n);
            Plin(l,l) = 0;
            Plin(l,i) = 1;
            Plin(i,i) = 0;
            Plin(i,l) = 1;
            U = Plin*U;
            L = L*Plin; % Plin^(-1) = Plin' = Plin pour une transposition
            P = Plin*P;
        end
        % Cas o� on intervertit 2 colonnes
        if(m ~= i)
            signatureQ = signatureQ + 1;
            Pcol = eye(n);
            Pcol(m,m) = 0;
            Pcol(i,m) = 1;
            Pcol(i,i) = 0;
            Pcol(m,i) = 1;
            U = U*Pcol;
            Q = Q*Pcol;
        end    
    end
    Li = eye(n);
    Li((i+1):n,i) = - U((i+1):n,i) / U(i,i);
    U = Li*U;
    L = L*(2*eye(n) - Li); %On utilise que L(i)^-1 = 2Ident - L(i)
end

L = P*L;

% Calcul du d�terminant de A
det = (-1)^(signatureQ+signatureP) * prod(diag(U));

if(det ~= 0)
       disp("La matrice est inversible : le produit est donc unique.");

end
end