%% Orthogonal Matching Pursuit

function [x_est] = OMP3(y_mat, eps, n, L, F_n, F_L)
    disp('OMP Initiated')
    K = 0;
    x_est = zeros(n, L);
    
    res_mat = y_mat;
%     basis_matrix = [];
    error = norm(res_mat);
%     x_k_hat = [];
%     support   = [];
    
    F1 = F_L';
    F2 = conj(F_n);
    
    while (error > eps)
        K = K+1;
        corr_mat = F_n * res_mat * transpose(F_L);
        [~, max_idx] = max(abs(corr_mat),[],'all');
        row_idx = 1+mod(max_idx-1, n);
        col_idx = 1+floor((max_idx-1)/n);
%         col = kron(F1(:, col_idx), F2(:, row_idx));
        b_k = corr_mat(row_idx, col_idx);
        x_est(row_idx, col_idx) = b_k;
%         x_k_hat = [x_k_hat; b_k];
%         basis_matrix = [basis_matrix col];
%         support = [support max_idx];
        res_mat = res_mat - (b_k * (F2(:, row_idx) * transpose(F1(:, col_idx))));
        error = norm(res_mat);
%         disp(error);
    end 
    
%     x(support) = x_k_hat;
    
    disp('Sparsity of OMP Result:')
    disp(length(find(x_est~=0)))
    disp('Iterations taken:')
    disp(K)
end