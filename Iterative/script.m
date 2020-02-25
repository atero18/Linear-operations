%{ 

EXPLICATIONS SUR LES FONCTIONS

creatAB cr�� la matrice A et le vecteur B demand� dans le TP (fonction de n)
partiesMat subdivise la matrice A en sa partie D,E,F
Les fonctions Gauss_Seidel, Jacobi et relaxation ne servent qu'� renvoyer
la matrice d'it�ration et M^-1

iterSolv s'occupe de toute la partie algorithmique

recherche_omega_opti recherche le meilleur omega en balayant l'intervalle
]0,2[
wOptiHerPos renvoi l'om�ga optimal pour une matrice hermitienne def.
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

% A est sym�trique dans R (donc Hermitienne) et toutes ses valeurs propres
% sont positives (eig(A)). Donc on peut calculer le omega optimal :
wOpti = wOptiHermPos(A);

% On peut aussi ne pas utiliser cette propri�t� et effectivement
% repr�senter le rayon spectral de la matrice en fonction de omega :
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

%Evolution du nombre d'it�rations en fonction du rayon spectral
plot(Data(3,:),Data(2,:));
xlabel("Rayon spectral");
ylabel("Nombre d'it�rations");

%% Comparaison du nombre d'it�rations pour les diff�rentes m�thodes

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
ylabel("Nombre d'it�rations");

%% Comparaison nombre it�rations empirique et minimal pour Jacobi

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
xlabel("rho (matrice it�ration)");
ylabel("Nombre d'it�rations");

%% Partie 2

B = [4 3 3; 3 4 3; 3 3 4] ./ 4;
disp("Valeurs propres de B");
eig(B)
% La matrice est sym�trique et ses v.p. sont strictements positives. Elle
% est donc d�finie positive.
% N�anmoins il n'y a pas convergence. En effet le rayon spectral de la
% matrice d'it�ration est sup�rieure � 1 :
Ite = Jacobi(B);
disp("Rayon spectral de la matrice d'it�ration");
max(abs(eig(Ite)))

C = [1 2 -2; 1 1 1; 2 2 1];
% On calcule le rayon spectral de la matrice d'it�ration
Ite = Jacobi(C);
max(abs(eig(Ite)))
disp("Le rayon spectral est inf�rieur � 1. Jacobi converge pour C");

Ite = Gauss_Seidel(C);
max(abs(eig(Ite)))
disp("Le rayon spectral est sup�rieur � 1. Gauss-Seidel diverge pour C");

D = [2 -1 1; 2 2 2;-1 -1 2];
Ite = Jacobi(D);
max(abs(eig(Ite)))
disp("Le rayon spectral est sup�rieur � 1. Jacobi diverge pour C");

Ite = Gauss_Seidel(D);
max(abs(eig(Ite)))
disp("Le rayon spectral est inf�rieur � 1. Gauss_Seidel converge pour C");    