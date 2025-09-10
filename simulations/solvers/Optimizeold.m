 function M = Optimizeold(Sigma_hat, gamma_j, alpha, ind, Llist)
        % Input:
        % Sigma_hat: The estimated covariance matrix (hat{Sigma})
        % gamma_j: The constraint parameter gamma_j
        % alpha: The threshold parameter for H_alpha
        % ind: the index matrix to define the groups
    
        % Dimensions
        [p, ~] = size(Sigma_hat);
    
        % Initialize the output matrix M
        M = zeros(p);
    
        % Number of elements in Llist
        TransL = numel(Llist);
    
        % Preallocate a cell array to store results from parfor
        M_temp = cell(TransL, 1);
    
        % Options for fmincon with parallel processing enabled
    options = optimoptions('fmincon', 'Algorithm', 'interior-point', ...
                           'Display', 'Iter', 'TolFun', 1e-9, 'TolCon', 1e-9, ...
                           'TolX', 1e-9, 'MaxIter', 30,'MaxFunctionEvaluations', 1e6);
    
        % Solve for each column of M (each l-th basis vector e_l) in parallel
        for i = 1:TransL
            l = Llist(i);
            disp(l);
    
            % Define the natural basis vector e_l
            e_l = zeros(p, 1);
            e_l(l) = 1;
    
            % Initial guess for m_opt
            m0 = zeros(p, 1);
    
            % Objective function for fmincon: (0.5 * m' Σ_hat * m)
            objective = @(m) 0.5 * m' * Sigma_hat * m;
    
            % Nonlinear constraint function
            nonlincon = @(m) group_norm_constraint(Sigma_hat * m - e_l, gamma_j, alpha, ind);
    
            % Solve the constrained optimization problem using fmincon
            m_opt = fmincon(objective, m0, [], [], [], [], [], [], nonlincon, options);
            disp(groupnorm(H_alpha(Sigma_hat * m_opt - e_l,alpha), inf, 2, ind));
            % Store the result in the temporary cell array
            M_temp{i} = m_opt;
        end
    
        % After parfor, assign the results back to M
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
