FID = fopen('adc_data.bin','r');
data = fread(FID,'int16');
fclose(FID);

n = 256;
L = 128;


output = data;
output_cmplx = zeros(256*128*8,4);

for ch=1:4
    for i = 1:128*256*8
        output_cmplx(i,ch) = complex(output(8*(i-1)+2*ch-1) , output(8*(i-1)+2*ch));             
    end
end

Rx_complex_signal = output_cmplx(:,1);

y = zeros(256*128, 1);

for i=1:256*128
    y(i) = Rx_complex_signal(i);
end

Y = reshape(y, n, L);

F_n = 1/sqrt(n) * dftmtx(n);
F_L = 1/sqrt(L) * dftmtx(L);

X_true = F_n * Y * transpose(F_L); 
x_true_vec = X_true(:);

figure(1);
surfc(abs(X_true), 'EdgeColor', 'none');
title("2D FFT of true ADC data");
xlabel("Doppler Axis");
ylabel("Range Axis");

% k = 23; % sparsity
threshold = 25000;

X_estimated = OMP3(Y, threshold, n, L, F_n, F_L);
x_est_vec = X_estimated(:);


figure(2);
surfc(abs(X_estimated), 'EdgeColor', 'none');
title("Estimated 2D FFT from 2D OMP algorithm");
xlabel("Doppler Axis");
ylabel("Range Axis");

% figure(2);
% subplot(3,1,1);
% stem(abs(X_true), 'og');
% title('X True');
% subplot(3,1,2);
% stem(abs(X_estimated), 'ob');
% title('X Estimated');
% subplot(3,1,3);
% stem(abs(X_true-X_estimated), 'or');
% title('Difference (True-Estimated)');