% Function to solve linear system using Jacobi method. It takes input A, b,
% n, initial estimate and tolerance. It returns the final estimate, no. of
% iterations for convergence and an array which contains error value
% obtained at the end of each iteration.
function[x, k, err_plt] = jacobi_method(A, b, n, x0, tolerance)
    x = x0;
    y = x0;
    err = norm(A*x-b, 2)/norm(A*x0-b, 2);
    k = 0;
    err_plt = [];
    while err >= tolerance
        for i=1:n
            t = 0;
            for j=1:n
                if i ~= j
                    t = t + A(i, j) * x(j);
                end
            end
            y(i) = (1/A(i, i)) *(-t + b(i));
        end
        err_plt = [err_plt norm(y-x, 2)/norm(x, 2)];
        x = y;
        k = k + 1;
        err = norm(A*x-b, 2)/norm(A*x0-b, 2);
    end
end
