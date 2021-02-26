clc;
clear all;
close all;

n = input('Enter size of matrix - ');
% Run for higher number of iterations for more accurate results
iter = input('Enter number of iterations - ');
% Generate a symmetric matrix
A = rand(n, n);
A = (A + transpose(A))/2;
% hessenberg is a function used to convert matrix to a hessenberg
% matrix
H = hessenberg(A, n);
for i=1:iter
    H0 = H;
    % givens function uses Givens rotation to find QR decomposition of a tridiagonal
    % matrix
    [Q, R] = givens(H0, n);
    H = R*Q;
end
% Eigen-values are the diagonal elements of Hessenberg matrix
eig_val = sort(diag(H));
eig_val_in_built = sort(eig(A));
% Sorted them in ascending order so that it becomes easier to find error in
% calculations
display(eig_val);
display(eig_val_in_built);
error = norm(eig_val - eig_val_in_built, 2);
fprintf('Error in calculations is %d\n\n', error);
fprintf('Per iteration complexity of the QR algorithm is O(n^4).\n');
fprintf('This is because there is one loop in the givens code which runs for n-1\n');
fprintf('times and in each iteration of it, there is matrix multiplication\n');
fprintf('whose complexity is O(n^3). Hence, total algorithm per iteration\n');
fprintf('complexity becomes O(n^4).\n');
