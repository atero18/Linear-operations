function [x,r,iter] = iterSolv(data,A,b,x0,eps,nbiter)

% Si 'data' est une fonction on suppose que c'est une méthode sans
% arguments (Jacobi, Gauss_Seidel..)
if isa(data,'function_handle')
    [matIter,invM] = data(A);

% Sinon c'est une cellule. Auquel cas soit c'est une fonction avec
% arguments, soit ce sont les matrices d'iteration et M^-1
else
    if isa(data{1}, 'function_handle')
        [matIter,invM] = data{1}(A,data{2});
    else
        matIter = data{1};
        invM = data{2};
    end
        
end

iter = 0;
x = x0;
r = norm(b - A*x);
while iter < nbiter && r > eps
    x = matIter * x + invM*b;
    r = norm(b - A*x);
    iter = iter + 1;
end
end

