% Helper function that returns the compressed image
function[J] = compression(aa, bb, x_R, x_G, x_B, mu_R, mu_G, mu_B, V_R, V_G, V_B, K)
    a = aa/8;
    b = bb/8;
    N = a * b;
    J = zeros(aa, bb, 3);
    x_tilde_R = zeros(64, N);
    x_tilde_G = zeros(64, N);
    x_tilde_B = zeros(64, N);
    for i=1:N
        x_tilde_R(:, i) = mu_R;
        x_tilde_G(:, i) = mu_G;
        x_tilde_B(:, i) = mu_B;
        for k=1:K
            x_tilde_R(:, i) = x_tilde_R(:, i) + dot(x_R(:, i) - mu_R, V_R(:, k)) * V_R(:, k);
            x_tilde_G(:, i) = x_tilde_G(:, i) + dot(x_G(:, i) - mu_G, V_G(:, k)) * V_G(:, k);
            x_tilde_B(:, i) = x_tilde_B(:, i) + dot(x_B(:, i) - mu_B, V_B(:, k)) * V_B(:, k);
        end
    end
    R_c = zeros(aa, bb);
    G_c = zeros(aa, bb);
    B_c = zeros(aa, bb);
    k = 1;
    for i=0:a-1
        for j=0:b-1
            y_r = x_tilde_R(:, k);
            y_r = reshape(y_r, [8, 8]);
            R_c(8*i+1:8*(i+1), 8*j+1:8*(j+1)) = y_r;
            y_g = x_tilde_G(:, k);
            y_g = reshape(y_g, [8, 8]);
            G_c(8*i+1:8*(i+1), 8*j+1:8*(j+1)) = y_g;
            y_b = x_tilde_B(:, k);
            y_b = reshape(y_b, [8, 8]);
            B_c(8*i+1:8*(i+1), 8*j+1:8*(j+1)) = y_b;
            k = k + 1;
        end
    end
    J(:, :, 1) = R_c;
    J(:, :, 2) = G_c;
    J(:, :, 3) = B_c;
end
