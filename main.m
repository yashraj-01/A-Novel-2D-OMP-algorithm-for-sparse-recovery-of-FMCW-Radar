clc;

n = 128;
L = 32;
k = 10; % sparsity
threshold = 1.0:0.01:7.5;

noise_var = 0.001;
noise = sqrt(noise_var) * randn(n, L);

p = randperm(n*L, k);
x_true = zeros(n*L, 1);
x_true(p) = (5*rand(k,1)-.5) + (5*rand(k,1)-.5)*1i;

x_mat = reshape(x_true, n, L);

x_mat = x_mat + noise;

F_n = 1/sqrt(n) * dftmtx(n);
F_L = 1/sqrt(L) * dftmtx(L);

y_mat = F_n' * x_mat * conj(F_L);

y = y_mat(:);

abs_diff = [];

for eps = threshold
   x_estimated = OMP3(y_mat, eps, n, L, F_n, F_L); 
   abs_diff = [abs_diff norm(abs(x_mat(:)-x_estimated(:)))];
end

figure(1);
plot(1./threshold, abs_diff);
xlabel("1/Threshold");
ylabel("Error");
title("Accuracy plot of optimized 2D OMP algorithm");
% figure;
% subplot(3,1,1);
% stem(abs(x_mat(:)), 'or');
% title('X True');
% subplot(3,1,2);
% stem(abs(x_estimated(:)), 'ob');
% title('X Estimated');
% subplot(3,1,3);
% stem(abs(x_mat(:)-x_estimated(:)), 'og');
% title('Difference (True-Estimated)');
