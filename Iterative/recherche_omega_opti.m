function [omegamin] = recherche_omega_opti(A)

n = size(A);
n = n(1);
[D,E] = partiesMat(A);

w = linspace(0.01,2,100);
rho = zeros(1,length(w));

for i = 1:length(w)
    rho(i) = max(abs(eig(eye(n) - inv(D ./ w(i) - E)*A)));
end

plot(w,rho);
xlabel("omega");
ylabel("rho");
[~,I] = min(rho);
omegamin = w(I(1));