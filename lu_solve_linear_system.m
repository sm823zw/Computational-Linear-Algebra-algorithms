function[x] = lu_solve_linear_system(L, U, b, n)
    % Forward Substitution
    z = zeros(n, 1);
    for i=1:n
        S = 0;
        for j=1:i-1
            S = S + L(i, j) * z(j);
        end
        z(i) = (1/L(i, i)) * (b(i) - S);
    end
    % Backward Substitution
    x = zeros(n, 1);
    for i=n:-1:1
        S = 0;
        for j=i+1:n
            S = S + U(i, j) * x(j);
        end
        x(i) = (1/U(i, i)) * (z(i) - S);
    end
end
