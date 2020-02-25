function [ray,n] = ray_spec(A,N)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin == 1
    N = 256;
end

n = 1;

%Calcul de A^N
while n < N
    A = A * A;
    n = 2 * n;
end

%Calcul du rayon spectral
norme1 = max(sum(abs(A)));
ray = norme1^(1/n);
end

