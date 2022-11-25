%% Orthogonal Matching Pursuit

function [x] = OMP1(y, eps, n, L, F_n, F_L)
    disp('OMP Initiated')
    K = 0;
    x = zeros(n*L, 1);
    
    res = y;
    basis_matrix = [];
    error = norm(res);
    x_k_hat = [];
    support   = [];
    
    F1 = F_L';
    F2 = conj(F_n);
    
    while error > eps
        K = K+1;
        res_mat = reshape(res,n,L);
        corr_mat = F_n * res_mat * transpose(F_L);
        [~, max_idx] = max(abs(corr_mat),[],'all');
        row_idx = 1+mod(max_idx-1, n);
        col_idx = 1+floor((max_idx-1)/n);
        col = kron(F1(:, col_idx), F2(:, row_idx));
        x_k_hat = [x_k_hat; corr_mat(row_idx, col_idx)];
        basis_matrix = [basis_matrix col];
        support = [support max_idx];
        res = y-basis_matrix*x_k_hat;
        error = norm(res);
    end 
    
    x(support) = x_k_hat;
    
    disp('Sparsity of OMP Result:')
    disp(length(find(x~=0)))
    disp('Iterations taken:')
    disp(K)
end