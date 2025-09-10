function M = OptimizeMnew(U,D,Wj, gamma_j, alpha, ind,Llist)
    % Input:
    % Sigma_hat: The estimated covariance matrix (hat{Sigma})
    % gamma_j: The constraint parameter gamma_j
    % alpha: The threshold parameter for H_alpha
    % ind: The index matrix to define the groups
    % Llist: List of indices for which to compute the optimization

    % Dimensions

    [n,p]=size(Wj);
    
    
    V=Wj'*U*diag((diag(D)).^(-0.5))/sqrt(n);
    

    M = zeros(p);

    % Number of elements in Llist
    TransL = numel(Llist);

    % Preallocate a cell array to store results
    M_temp = cell(TransL, 1);

    % Options for fmincon
    options = optimoptions('fmincon', 'Algorithm', 'interior-point', ...
        'Display', 'off', 'TolFun', 1e-6, 'TolCon', 1e-6, ...
        'TolX', 1e-6, 'MaxIter', 80, 'MaxFunctionEvaluations', 1e6,'SpecifyObjectiveGradient', true, ...
    'SpecifyConstraintGradient', false, ... 
    'CheckGradients', false);

    % Precompute V' * e_l for all l in Llist
    % Vt_e_l_all = zeros(n, TransL);
    % e_l_perp_all = zeros(p, TransL);
    % for i = 1:TransL
    %     l = Llist(i);
    %     e_l = zeros(p, 1);
    %     e_l(l) = 1;
    %     % Vt_e_l = V' * e_l;           % n x 1 vector
    %     % Vt_e_l_all(:, i) = Vt_e_l;

    %     % % Compute the orthogonal component e_l_perp = e_l - V * (V' * e_l)
    %     % e_l_perp = e_l - V * Vt_e_l; % p x 1 vector
    %     % e_l_perp_all(:, i) = e_l_perp;
    % end
    vecD=diag(D);
    % Solve for each column of M
    % if isempty(gcp('nocreate'))  % Check if a parallel pool is already running
    %     parpool('local', 10);     % Start the pool with 9 workers
    % end
    for i = 1:TransL
        l = Llist(i);
        disp(['Processing index: ', num2str(l)]);
        e_l = zeros(p, 1);
        e_l(l) = 1;

        % Extract V' * e_l and e_l_perp for the current l
        % Vt_e_l = Vt_e_l_all(:, i);
        % e_l_perp = e_l_perp_all(: , i);

        % Initial guess for u_opt (n-dimensional)
        u0 = ones(n, 1)/sqrt(n);

        % Objective function in terms of u (n-dimensional)
        objective = @(u) deal( 0.5 * sum(vecD .* (u.^2)), vecD .* u);

        % Nonlinear constraint function in terms of u
        % Compute s = D .* u - Vt_e_l (n-dimensional)
        % Compute m_effective = V * s + e_l_perp (p-dimensional)
        nonlincon = @(u) group_norm_constraint(V * (vecD .* u) - e_l, gamma_j, alpha, ind);

        % Solve the constrained optimization problem using fmincon
        [u_opt, ~, exitflag] = fmincon(objective, u0, [], [], [], [], [], [], nonlincon, options);

        % Check if the optimization was successful
        if exitflag <= 0
            warning(['Optimization did not converge for index ', num2str(l), '.']);
        end

        % Recover m_opt from u_opt
        % u_full = [u_opt; zeros(p - n, 1)];
        m_opt = V * u_opt;

        % Store the result in the temporary cell array
        M_temp{i} = m_opt;

        % disp(group_norm_constraint(V * (vecD .* u_opt) - e_l, gamma_j, alpha, ind));
        % disp(groupnorm(H_alpha(Sigma_hat * m_opt - e_l,alpha), inf, 2, ind));
    end
    % delete(gcp('nocreate'));
    % Assign the results back to M
    for i = 1:TransL
        l = Llist(i);
        M(:, l) = M_temp{i};
    end
end

function H_alpha_x = H_alpha(x, alpha)
    % Soft-thresholding function H_alpha
    H_alpha_x = sign(x) .* max(abs(x) - alpha, 0);
end

function [c, ceq] = group_norm_constraint(m_effective, gamma_j, alpha, ind)
    % Apply the soft-thresholding function H_alpha
    H_alpha_val = H_alpha(m_effective, alpha);

    % Compute the group norm
    norm_value = compute_group_norm(H_alpha_val, inf, 2, ind);

    % Inequality constraint: group norm <= gamma_j
    c = norm_value - gamma_j;
 
    % No equality constraints
    ceq = [];
end

function norm_value = compute_group_norm(gamma, q1, q2, ind)
    % Compute the q1, q2 group norm function
    num_groups = size(ind, 2);
    group_norm_values = zeros(1, num_groups);

    % Compute the q2 norm for each group
    for j = 1:num_groups
        start_idx = ind(1, j);
        end_idx = ind(2, j);

        % Extract the j-th group from gamma
        gamma_j = gamma(start_idx:end_idx);

        % Compute the q2-norm of the group
        group_norm_values(j) = norm(gamma_j, q2);
    end

    % Compute the q1-norm across all group norms
    norm_value = norm(group_norm_values, q1);
end
