clc;
clear all;
close all;

n = input('Enter size of the matrix - ');
iter = input('Enter number of iterations - ');
% Run for higher number of iterations for more accurate results
A = randn(n, n);
B = transpose(A) * A;
H = hessenberg(B, n);
for i=1:iter
    H0 = H;
    % givens function uses Givens rotation to find QR decomposition of a tridiagonal
    % matrix
    [Q, R] = givens(H0, n);
    H = R*Q;
end
% Singular values are square root of the eigen-values which are in turn 
% the diagonal elements of Hessenberg matrix
sv = sort(sqrt(diag(H)));
sv_in_built = sort(svd(A));
% Sorted them in ascending order so that it becomes easier to find error in
% calculations
display(sv);
display(sv_in_built);
error = norm(sv - sv_in_built, 2);
fprintf('Error in calculations is %d\n', error);
