function generate(n,p,q,nu,Treptition)
    dim=(p-1)*(q+1); % dimension in each nodewise regression


%% true coefficient represented as a p*p*(q+1) array
tB = tensor(zeros(1,p*p*(q+1)), [p p q+1]); 

listi=2:3;
listj=[1];
for i=listi
	for j=listj
		tB(j,j+1,i) = nu;
        tB(j+1,j,i) = tB(j,j+1,i);
    end
end
path=append("datarepo/n=",string(n),"p=",string(p),"q=",string(q),"nu=",string(nu),'/');
mkdir(path);

beta = [];

% Loop over each j to construct beta_j
for j = 1:p
    beta_j = [];  % Initialize empty beta_j vector
    
    % Loop over each k and h
    for h = 1:(q+1)
        for k = 1:p
            if k ~= j  % Exclude cases where j == k
                % Append the value -tB(j, k, h) to beta_j
                element=tB(j, k, h);
                beta_j = [beta_j; -tB(j, k, h)];
            end
        end
    end
    
    % Concatenate beta_j to the full beta vector
    beta = [beta; beta_j];
end
csvwrite(append(path,'beta.csv'),beta)

% while reptition <= 1000
%     N = randi([0 1], n, q);
    
%     X = zeros(n, p);  % Initialize X
%     a = 0;  % Initialize counter for successfully generated X rows

%     for i = 1:n
%         omega = ttv(tB, iU(i,:)', 3);
%         omega = double(omega);
%         omega = omega + 1.0 * diag(ones(1, p));  % Add identity matrix to omega
%         sigma = inv(omega);

%         % Initialize positive definiteness flag
%         isPositiveDefinite = true;
        
%         % Check if omega is positive definite
%         try
%             chol(omega);  % Try Cholesky decomposition
%         catch
%             isPositiveDefinite = false;  % If it fails, mark as not positive definite
%         end

%         if isPositiveDefinite
%             % Generate data
%             X(i, :) = mvnrnd(zeros(1, p), sigma, 1);
%             a = a + 1;  % Successfully generated one row of X, increment counter
%         else
%             break;  % If omega is not positive definite, break and retry
%         end
%     end
    
%     % If all n rows of X are successfully generated
%     if a == n
%         % Create directory and save data
%         path = append('datarepo/rep', string(reptition), '/');
%         mkdir(path);
%         csvwrite(append(path, 'iU.csv'), iU);
%         csvwrite(append(path, 'Z0.csv'), X);
%         reptition = reptition + 1;  % Increment reptition counter
%     end
% end

for reptition=1:Treptition
    path0=append(path,'/rep',string(reptition),'/');
    mkdir(path0);
    N = randi([0 1],n,q);
    U = (N-0.5)/0.5;  %对列做归一化
    iU = [ones(n,1), U];
    csvwrite(append(path0,'iU.csv'),iU);
    X = zeros(n, p);

    for i=1:n
        omega = ttv(tB,iU(i,:)', 3);
        omega = double(omega);
        omega = omega + 1.0*diag(ones(1,p));
        sigma = inv(omega);
        X(i,:) = mvnrnd(zeros(1,p),sigma,1);
    end
csvwrite(append(path0,'Z0.csv'),X);
end

end