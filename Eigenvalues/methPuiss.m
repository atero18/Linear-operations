function [ray,vecP,lambdak,Normex] = methPuiss(A,x0,N,tol)


% On peut choisir de ne pas mettre en place uniquement une N-itération sans
% prendre en compte la tolérance (dans le cas par exemple où x0 est
% orthogonal à E(lambda)

if(nargin == 3)
    tol = -2;
end

dim = size(A);
dim = dim(1);


%Si le nombre d'itérations est fini on créé directement
% les tableaux pour optimiser la mémoire

if N < Inf
    K = N + 1;
else
    K = 2;
end

x = zeros(dim,K);
q = zeros(dim,K);
lambdak = zeros(dim,K);
Normex = zeros(1,dim);
x(:,1) = x0;

Normex(1) = norm(x0);
q(:,1) = x0 / Normex(1);
% Lambda(0) n'existe pas. On le considère comme x(1) / q(1)
lambdak(:,1) = norm(1);

i = 1;
delta = tol + 1;
while i <= N && delta >= tol
    i = i+1;
    x(:,i) = A * q(:,i-1);
    Normex(i) = norm(x(:,i));
    
    lambdak(:,i) = x(:,i) ./ q(:,i-1);
    %Cas où une composante du vecteur q(i-1) serait nulle.
    for j = 1:dim
        if(norm(q(j,i-1)) < eps)
            lambdak(j,i) = sign(x(j,i))*10^8;
        end
    end
    q(:,i) = x(:,i) / Normex(i);
    
    delta = norm(lambdak(:,i) - lambdak(:,i-1),Inf);
    
end

%Calcul de la valeur propre et du vecteur propre
valP = min(lambdak(:,i));
ray = norm(valP);
vecP = (conj(valP) / ray)^i * q(:,i);
end

