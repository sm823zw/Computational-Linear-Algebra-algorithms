% LU Decomposition using partial pivoting
function[P, L, U] = lu_decomposition(A, n)
    I = eye(n);
    U = A;
    L = I;
    P = I;
    for i = 1:n-1
        piv = max(abs(U(i:n, i)));
        for j = i:n
            if abs(U(j, i)) == piv
                piv_row_index = j;
                break;
            end
        end
        temp = U(piv_row_index, i:n);
        U(piv_row_index, i:n) = U(i, i:n);
        U(i, i:n) = temp;
        temp = L(piv_row_index, 1:i-1);
        L(piv_row_index, 1:i-1) = L(i, 1:i-1);
        L(i, 1:i-1) = temp;
        temp = P(piv_row_index, :);
        P(piv_row_index, :) = P(i, :);
        P(i, :) = temp;
        for j = i+1:n
            alpha = U(j, i) / U(i, i);
            L(j, i) = alpha;
            U(j, :) = U(j, :) - alpha * U(i, :);
        end
    end
end
