% QR Decomposition using Gram-Schmidt procedure
function[Q, R] = gs(A, n)
    Q = zeros(n);
    for i = 1:n
        a = A(:, i);
        q = a; 
        if i ~= 1
            for j = 1:i
                q = q - (dot(a, Q(:, j)))*Q(:, j);
            end
        end
        q = q / norm(q, 2); % Normalizing the vector
        Q(:, i) = q;
    end
    R = transpose(Q) * A;
end
