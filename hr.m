% QR Decomposition using Householder reflections
function[Q, R] = hr(A, n)
    I = eye(n);
    Q = eye(n);
    W = A;
    for i = 1:n
        a = W(i:n, i);
        norm_a = sqrt(sum(a.^2));
        if sign(a(1)) >= 0
            v = a + norm_a * I(i:n, i);
        else
            v = a - norm_a * I(i:n, i);
        end
        h = eye(n - i + 1) - 2 * (v * transpose(v)) / (transpose(v) * v);
        H = I;
        H(i:n, i:n) = h;
        Q = Q * H;
        W = H * W;
    end
    R = W;
end
