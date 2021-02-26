% QR Decomposition using rotation matrix tailor-made for tridiagonal
% matrices
function[Q, R] = givens(A, n)
    I = eye(n);
    Q = eye(n);
    W = A;
    for i = 1:n-1
        j = i+1;
        v = W(j-1:j, i);        
        c = sqrt(v(1)^2 + v(2)^2);
        cos_t = v(1)/c;
        sin_t = v(2)/c;
        r = [cos_t -sin_t; sin_t cos_t];
        H = I;
        H(j-1:j, j-1:j) = r;
        W = transpose(H) * W;
        Q = Q * H;
    end
    R = W;
end
