clc;
clear all;
close all;
N = [10, 50, 100]; % Dimensions of matrix
for q=1:3
    n = N(q);
    e = ones(n,1);
    A = spdiags([-e 2*e -e], -1:1, n, n);
    A = full(A);
    b = rand(n,1);
    epsilon = [0.1, 0.01, 0.001, 0.0001, 0.00001];
    error = [];
    for ep=1:5
        % Jacobi method
        x0 = zeros(n, 1);
        [x, k] = jacobi_method(A, b, n, x0, epsilon(ep));
        e = norm(x - A\b, 2);
        error = [error e];
        fprintf('Error(with n = %d and epsilon = %d) = %f\n', n, epsilon(ep), e);
        fprintf('No. of iterations(with n = %d and epsilon = %d) = %d\n', n, epsilon(ep), k);
    end
    figure;
    plot(categorical(epsilon), error);
    title(['Error plot for n = ',num2str(n),'']);
    xlabel('Tolerance');
    ylabel('Error');
    fprintf('\n');
end
fprintf('As seen from the graph, for a fixed value of n, the error reduces if we reduce the tolerance.\n')
fprintf('The error reduces by the same factor we reduce the tolerance by.\n')
fprintf('Since, we reduce the tolerance by a factor of 10, the error also reduces by a factor of almost 10.\n')
fprintf('The error is larger as we increase the size of matrices from 10 to 100.\n');
