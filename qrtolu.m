clc;
clear all;
close all;
%%
n = input('Enter dimension of the matrix - ');
% A = randn(n, n) + randn(n, n)*1i;
A = randn(n, n);
I = eye(n);
[Q, R] = hr(A, n);

% Inverse of the matrix calculations
% inv_A = R\transpose(Q);
Z = transpose(Q);
% Backward Substitution
inv_A = zeros(n, n);
for i=n:-1:1
    S = 0;
    for j=i+1:n
        S = S + R(i, j) * inv_A(j, :);
    end
    inv_A(i, :) = (1/R(i, i)) * (Z(i, :) - S);
end
fprintf('Error in inverse calculation is %d\n', norm(A * inv_A - I, 'fro'));

% Determinant of the matrix calculations

% Since, determinant of Q is 1 or -1 and we don't know the sign,
% we can calculate the magnitude of the determinant of A.
mag_det_Q = 1;
det_R = 1;
for i=1:n
    det_R = det_R * R(i, i);
end

mag_det_A = mag_det_Q * det_R;
fprintf('Magnitude of determinant of the matrix is %d\n', abs(mag_det_A));

% Magnitude of determinant of A
det_A = det(transpose(Q)) * det_R;
fprintf('Determinant of the matrix is %d\n', det_A);
fprintf('Determinant of the matrix by in-built method is %d\n', det(A));

%%
n = input('Enter the dimension of the matrix - ');
A = randn(n, n);
[Q, R] = hr(A, n);

% Given QR, finding LU without constructing matrix A.
[PQ, LQ, UQ]= lu_decomposition(Q, n);
% Permutation matrix will be the same.
P = PQ;
% LQ is the lower triangular matrix
L = LQ;
% UQ and R are both upper triangular matrix. Hence, their multiplication
% will also be an upper triangular matrix.
U = UQ*R;

% Error calculations
E = norm(A - transpose(P)*L*U, 'fro');
fprintf('Error in calculations is %d\n', E);
