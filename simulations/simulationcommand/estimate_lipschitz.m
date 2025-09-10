function Lip = estimate_lipschitz(A, num_iter)
    % 估算矩阵A的最大奇异值（即Lipschitz常数）
    if nargin < 2
        num_iter = 20;
    end
    [n, ~] = size(A);
    x = randn(n,1);
    x = x/norm(x);
    for i=1:num_iter
        x = A * (A' * x);
        x = x / norm(x);
    end
    Lip = norm(A * (A' * x)) / norm(x);  % 这是最大奇异值
end