% function performs Cholesky decomposition and returns a lower triangular
% matrix
function[L] = cholesky(A, n)
    L = zeros(n, n);
    for i=1:n
        if i==1
            L(i, i) = sqrt(A(i, i));
        else
            L(i, 1) = A(i, 1)/L(1, 1);
        end
    end
    for i=2:n
        L(i, i) = A(i, i);
        for k=1:i-1
            L(i, i) = L(i, i) - L(i, k)^2;
        end
        L(i, i) = sqrt(L(i, i));
        if i~=n
            for j=i+1:n
                L(j, i) = A(j, i);
                for k=1:i-1
                    L(j, i) = L(j, i) - L(i, k) * L(j, k);
                end
                L(j, i) = L(j, i)/L(i, i);
            end
        end
    end
end
