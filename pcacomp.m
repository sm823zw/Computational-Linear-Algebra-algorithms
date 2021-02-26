clc;
clear all;
close all;

% Separate the three channels
I = im2double(imread('watch.bmp'));
R = I(:, :, 1);
G = I(:, :, 2);
B = I(:, :, 3);
[aa, bb] = size(R);
a = aa/8;
b = bb/8;
N = a * b;

% Calculate the mean
x_R = zeros(64, N);
x_G = zeros(64, N);
x_B = zeros(64, N);
mu_R = zeros(64, 1);
mu_G= zeros(64, 1);
mu_B = zeros(64, 1);
k = 1;
for i=0:a-1
    for j=0:b-1
        x_r = R(8*i+1:8*(i+1), 8*j+1:8*(j+1));
        x_r = reshape(x_r, [64, 1]);
        mu_R = mu_R + x_r;
        x_R(:, k) = x_r;
        x_g = G(8*i+1:8*(i+1), 8*j+1:8*(j+1));
        x_g = reshape(x_g, [64, 1]);
        mu_G = mu_G + x_g;
        x_G(:, k) = x_g;
        x_b = B(8*i+1:8*(i+1), 8*j+1:8*(j+1));
        x_b = reshape(x_b, [64, 1]);
        mu_B = mu_B + x_b;
        x_B(:, k) = x_b;
        k = k + 1;
    end
end
mu_R = mu_R/N;
mu_G = mu_G/N;
mu_B = mu_B/N;

% Calculate covariance matrix
sigma_R = zeros(64, 64);
sigma_G = zeros(64, 64);
sigma_B = zeros(64, 64);
for i=1:N
    sigma_R = sigma_R + (x_R(:, i) - mu_R) * transpose(x_R(:, i) - mu_R);
    sigma_G = sigma_G + (x_G(:, i) - mu_G) * transpose(x_G(:, i) - mu_G);
    sigma_B = sigma_B + (x_B(:, i) - mu_B) * transpose(x_B(:, i) - mu_B);
end
sigma_R = sigma_R/N;
sigma_G = sigma_G/N;
sigma_B = sigma_B/N;

% Get the eigenvalues and eigenvectors
[V_R,eig_R] = eig(sigma_R);
[V_G,eig_G] = eig(sigma_G);
[V_B,eig_B] = eig(sigma_B);
eig_R = diag(eig_R);
eig_G = diag(eig_G);
eig_B = diag(eig_B);
[eig_R,indices_R] = sort(eig_R,'descend');
[eig_G,indices_G] = sort(eig_G,'descend');
[eig_B,indices_B] = sort(eig_B,'descend');
V_R = V_R(:, indices_R);
V_G = V_G(:, indices_G);
V_B = V_B(:, indices_B);

%%
% Compression
K = [1, 5, 10, 20, 64];
for i=1:5
    % compression is a helper function that returns the compressed image
    J = compression(aa, bb, x_R, x_G, x_B, mu_R, mu_G, mu_B, V_R, V_G, V_B, K(i));
    figure; imshow(J); title(['Compressed Image, K= ', num2str(K(i)), '']);
end
figure; imshow(I); title('Original Image');

%%
% Error plot
error = zeros(64, 1);
for i=1:64
    % compression is a helper function that returns the compressed image
    J = compression(aa, bb, x_R, x_G, x_B, mu_R, mu_G, mu_B, V_R, V_G, V_B, i);
    error(i) = sqrt(sum(sum(sum((I-J).*(I-J)))));
end
figure;
K = 1:64;
plot(K, error);
title('Error vs Dimension of Subspace');
xlabel('Dimension of Subspace');
ylabel('Error');

% Comments
fprintf('The error decreases monotonically as we increase the dimension of subspace.\n');
fprintf('The decrease in error is drastic as we increase the dimension of the subspace\n');
fprintf('from 1 to around 10. After that error decreases rather gradually.\n');
