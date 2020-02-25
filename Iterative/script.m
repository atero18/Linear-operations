%{ 

EXPLICATIONS SUR LES FONCTIONS

creatAB créé la matrice A et le vecteur B demandé dans le TP (fonction de n)
partiesMat subdivise la matrice A en sa partie D,E,F
Les fonctions Gauss_Seidel, Jacobi et relaxation ne servent qu'à renvoyer
la matrice d'itération et M^-1

iterSolv s'occupe de toute la partie algorithmique

recherche_omega_opti recherche le meilleur omega en balayant l'intervalle
]0,2[
wOptiHerPos renvoi l'oméga optimal pour une matrice hermitienne def.
positive

 
%}
 

%% Question 2
[A,B] = creatAB(10);
[~,r1] = iterSolv(@Jacobi,A,B,B,0,100);
[~,r2] = iterSolv(@Gauss_Seidel,A,B,B,eps,100);
[~,r3] = iterSolv({@relaxation,3/2},A,B,B,eps,100);

r1,r2,r3

%% Question 3

[A,B] = creatAB(20);

% A est symétrique dans R (donc Hermitienne) et toutes ses valeurs propres
% sont positives (eig(A)). Donc on peut calculer le omega optimal :
wOpti = wOptiHermPos(A);

% On peut aussi ne pas utiliser cette propriété et effectivement
% représenter le rayon spectral de la matrice en fonction de omega :
wGraph = recherche_omega_opti(A)
wOpti, wGraph
%% Question 4

i = 1;

% n, nbiterations, rho
Data = [0 ; 0 ; 0];

for n = 4:40
    Data(1,i) = n;
    [A,B] = creatAB(n);
    
    [w,Data(3,i)] = wOptiHermPos(A); 
    [x,r,iterJ] = iterSolv({@relaxation,w},A,B,B,10^-12,Inf);
    Data(2,i) = iterJ;
      
    i = i + 1;
end
clear x r iter n i w;
clear A B;

%Evolution du nombre d'itérations en fonction du rayon spectral
plot(Data(3,:),Data(2,:));
xlabel("Rayon spectral");
ylabel("Nombre d'itérations");

%% Comparaison du nombre d'itérations pour les différentes méthodes

nbiterJac = zeros(1,37);
nbiterGauss = zeros(1,37);
for n = 4:40
    [A,B] = creatAB(n);
    [~,~,iterJ] = iterSolv(@Jacobi,A,B,B,10^-12,Inf);
    [~,~,iterG] = iterSolv(@Gauss_Seidel,A,B,B,10^-12,Inf);
    nbiterJac(n-3) = iterJ;
    nbiterGauss(n-3) = iterG;
end
    clear iterJ iterG n;

plot(4:40,nbiterJac,4:40,nbiterGauss,4:40,Data(2,:));
legend('Jacobi','Gauss Seidel','Relaxation');
xlabel("Taille n de la matrice");
ylabel("Nombre d'itérations");

%% Comparaison nombre itérations empirique et minimal pour Jacobi

kMin = zeros(1,37);
raySpec = zeros(1,37);
for n = 4:40
    [A,B] = creatAB(n);
    [Miter,invM] = Jacobi(A);
    x1 = Miter * B + invM*B;
    raySpec(n-3) = abs(max(eig(Miter)));
    kMin(n-3) = 1 + (log(10^-12) - log(norm(x1 - B))) / (log(raySpec(n-3)));
end
clear Miter invM n x1

plot(raySpec,nbiterJac,raySpec,kMin);
legend("nombre empirique","nombre minimal");
xlabel("rho (matrice itération)");
ylabel("Nombre d'itérations");

%% Partie 2

B = [4 3 3; 3 4 3; 3 3 4] ./ 4;
disp("Valeurs propres de B");
eig(B)
% La matrice est symétrique et ses v.p. sont strictements positives. Elle
% est donc définie positive.
% Néanmoins il n'y a pas convergence. En effet le rayon spectral de la
% matrice d'itération est supérieure à 1 :
Ite = Jacobi(B);
disp("Rayon spectral de la matrice d'itération");
max(abs(eig(Ite)))

C = [1 2 -2; 1 1 1; 2 2 1];
% On calcule le rayon spectral de la matrice d'itération
Ite = Jacobi(C);
max(abs(eig(Ite)))
disp("Le rayon spectral est inférieur à 1. Jacobi converge pour C");

Ite = Gauss_Seidel(C);
max(abs(eig(Ite)))
disp("Le rayon spectral est supérieur à 1. Gauss-Seidel diverge pour C");

D = [2 -1 1; 2 2 2;-1 -1 2];
Ite = Jacobi(D);
max(abs(eig(Ite)))
disp("Le rayon spectral est supérieur à 1. Jacobi diverge pour C");

Ite = Gauss_Seidel(D);
max(abs(eig(Ite)))
disp("Le rayon spectral est inférieur à 1. Gauss_Seidel converge pour C");    