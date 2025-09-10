function norm_val = groupnorm(gamma, q1, q2, ind)
    num_groups = size(ind, 2); % Number of groups
    group_norms = zeros(1, num_groups); % Store norms for each group

    % Compute the q2 norm for each group
    for j = 1:num_groups
        start_idx = ind(1, j);
        end_idx = ind(2, j);

        % Extract the j-th group from gamma
        gamma_j = gamma(start_idx:end_idx);

        % Compute the q2-norm of the group
        group_norms(j) = norm(gamma_j, q2);
    end

    % Compute the q1-norm across all group norms
    norm_val = norm(group_norms, q1);
end