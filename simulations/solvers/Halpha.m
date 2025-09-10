function H_alpha_x = Halpha(x, alpha)
    H_alpha_x = sign(x) .* max(abs(x) - alpha, 0);
end