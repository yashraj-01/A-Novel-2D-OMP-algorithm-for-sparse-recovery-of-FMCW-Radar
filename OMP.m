%% Orthogonal Matching Pursuit

function [x] = OMP(y, eps, n, L, F_n, F_L)
    disp('OMP Initiated')
    K = 0;
    x = zeros(n*L, 1);
    
    res = y;
    basis_matrix = [];
    error = norm(res);
    x_k_hat = [];
    support   = [];
    
    y_mat = reshape(y,n,L);
    lookup_matrix = F_n * y_mat * transpose(F_L);
    lookup_vector = lookup_matrix(:);
    lookup = abs(lookup_vector);
    F1 = F_L';
    F2 = conj(F_n);
    
    while error > eps
        K = K+1;
        [~, b] = max(lookup);
        col = kron(F1(:, 1+floor((b-1)/n)), F2(:, 1+mod(b-1, n)));
        basis_matrix = [basis_matrix col];
        support = [support b];
        x_k_hat = [x_k_hat; lookup_vector(b)];
        lookup(b) = -1;
        res = y-basis_matrix*x_k_hat;
        error = norm(res);
    end 
    
    x(support) = x_k_hat;
    
    disp('Sparsity of OMP Result:')
    disp(length(find(x~=0)))
    disp('Iterations taken:')
    disp(K)
end