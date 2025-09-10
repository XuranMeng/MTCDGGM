
addpath(genpath(pwd));
savepath


Treptition=200;

startrep = 1;

p = 70;
q = 120;
% if n == 200
%     startrep = 1;
% else
%     startrep = 1;
% end

nu=0.3;
dim=(p-1)*(q+1);


parmain.maxiter = 30;       % maximum outer iterations
parmain.stoptol = 1e-6;     % stopping tolerance
parmain.printyes = 1;       % print progress
parmain.Lip = 1;            % Lipschitz constant (for example)
parmain.scale = 1;        % enable adaptive rescaling
parmain.dscale = ones(q,1); % no additional column scaling
parmain.Ascaleyes = 0;       % disable column scaling (set to 1 if desired)
parmain.bscale = 1;         % scaling factors (set to 1 initially)
parmain.cscale = 1;
parmain.sigma = 1;          % initial penalty parameter sigma
parmain.orgojbconst = 0;    % any constant to add to the objective

theta0 = zeros(q,1);  % primal variable x in R^n

y0 = zeros(q,1);  % dual variable y in R^n


if n==400
    pools=12;
elseif n==200
    pools = 13;
else
    pools=10;
end

parpool('local', pools) 
parfor reptition=startrep:Treptition
    disp([n p q reptition]);
    path=append("datarepounknown/n=",string(n),"p=",string(p),"q=",string(q),"nu=",string(nu),'/rep',string(reptition),'/');
    Z0=csvread(append(path,'Z0.csv'));
    iU=csvread(append(path,'iU.csv'));
    U=iU(:,2:end);
    A=U;
    hatbGamma=[];
    for j=1:p
        y=Z0(:,j);
        lambda_grid = logspace(-2, 1, 10);   % 从 0.01 到 10 之间取 10 个点

        K = 5;  
        cv_error = zeros(length(lambda_grid), 1);

        cv = cvpartition(n, 'KFold', K); 

        for i = 1:length(lambda_grid)
            lambda = lambda_grid(i);
            err = 0;
            for k = 1:K
                train_idx = training(cv, k);
                test_idx  = test(cv, k);

                A_train = A(train_idx, :);
                y_train = y(train_idx);

                A_test = A(test_idx, :);
                y_test = y(test_idx);

                n0=size(A_train,1);
                xi0 = zeros(n0,1); % dual multiplier ξ in R^m
                [~, ~, ~, hatgamma] = Classic_Lasso_SSNAL(A_train, y_train, q, sqrt(n0)* lambda, parmain, y0, xi0, theta0);


                y_pred = A_test * hatgamma;

                err = err + mean((y_test - y_pred).^2);
            end
            cv_error(i) = err / K; 
        end


[~, best_idx] = min(cv_error);
best_lambda = lambda_grid(best_idx);

xi0 = zeros(n,1);
[~, ~, ~, hatgamma_final] = Classic_Lasso_SSNAL(A, y, q, best_lambda*sqrt(n), parmain, y0, xi0, theta0);
hatbGamma=[hatbGamma,hatgamma_final];
end
csvwrite(append(path,'hatbGamma.csv'),hatbGamma);
end


delete(gcp('nocreate'));
