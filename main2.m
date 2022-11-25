[y_mat, n, L] = readBinFile('combined_if_signal.bin');
[x_true,~,~] = readBinFile('2D_fft_output.bin');
x_true_vec = x_true(:);

F_n = 1/sqrt(n) * dftmtx(n);
F_L = 1/sqrt(L) * dftmtx(L);

k = 1; % sparsity
threshold = 1e-1;

x_estimated = OMP2(y_mat, threshold, n, L, F_n, F_L, k);
x_est_vec = x_estimated(:);

figure;
subplot(3,1,1);
stem(abs(x_true_vec), 'or');
title('X True');
subplot(3,1,2);
stem(abs(x_est_vec), 'ob');
title('X Estimated');
subplot(3,1,3);
stem(abs(x_true_vec-x_est_vec), 'og');
title('Difference (True-Estimated)');