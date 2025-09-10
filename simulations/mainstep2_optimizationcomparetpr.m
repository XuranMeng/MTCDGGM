mex mexfun/mexbwsolve.c
mex mexfun/mexfwsolve.c
mex mexfun/mexfz.c
mex mexfun/mexMatvec.c
mex mexfun/mexProjL2.c
mex mexfun/mexProxL2.c
mex mexfun/mextriang.c
mex mexfun/mexMatvecmultiple.c
mex mexfun/mexmatrixmul.c
addpath(genpath(pwd));
savepath
p = 70; % dimension of response
q = 120; % dimension of covariates
nu=0.3;
dim=(p-1)*(q+1);
Treptition=1;
% if n==400
%     Treptition=1;
% elseif n==200
%     Treptition=1;
% else
%     Treptition=1;
% end
grpsize = p-1; 
clear ind;
for i = 1:(q+1)
    if i == 1
        ind(1,1) = 1; ind(2,1) = grpsize; ind(3,1) = 0;
    else
        ind(1,i) = ind(2,i-1) + 1;
        ind(2,i) = i*grpsize;
        ind(3,i) = 1;
    end
end
ind


path=append("datarepo/n=",string(n),"p=",string(p),"q=",string(q),"nu=",string(nu),'/');
    

% beta_j=beta((1:(p-1)*(q+1))+(j-1)*(p-1)*(q+1),:);
% hatbeta_j=hatbeta((1:(p-1)*(q+1))+(j-1)*(p-1)*(q+1),:);


% for j=1:2
%     beta_j=beta((1:(p-1)*(q+1))+(j-1)*(p-1)*(q+1),:);
%     Llist=find(beta_j~=0);
%     Xj=interM;
%     Xj(:,((0:q)*p+j)) = [];
%     y1=Z0(:,j);
%     b=y1-mean(y1);
%     M = OptimizeMparallel(Sigma_hat, 0.08, 0.04, ind,Llist);
% end

if isempty(gcp('nocreate'))  % Check if a parallel pool is already running
    parpool('local', 10);     % Start the pool with 9 workers
end
parfor reptition=1:200
    disp(reptition);
    path=append("datarepo/n=",string(n),"p=",string(p),"q=",string(q),"nu=",string(nu),'/rep',string(reptition),'/');
    Z0=csvread(append(path,'Z0.csv'));
    iU=csvread(append(path,'iU.csv'));
    if n==800
        beta=csvread(append(path,'hatbeta_mullbde0.3alpha2.csv'));
    else
        beta=csvread(append(path,'hatbeta_mulcv.csv'));
    end
    
    interM = [];
    for l=1:(q+1)
        interM=[interM Z0.*iU(:,l)]; %these are the interaction terms between Z0 and iU
    end

    % if ~exist(path, 'dir')
    %     mkdir(path); % 如果文件夹不存在则创建
    % end
    for j=1
        beta_j=beta((1:(p-1)*(q+1))+(j-1)*(p-1)*(q+1),:);
        Llist=find(beta_j~=0);
        Wj=interM;
        Wj(:,((0:q)*p+j)) = [];
        y1=Z0(:,j);
        b=y1-mean(y1);
        % Sigma_hat=mexmatrixmul(Wj',Wj)/n;
        Xi=mexmatrixmul(Wj,Wj')/n;
        [U,D]=eig(Xi);
   

        M = OptimizeMnew(U, D, Wj, 2/sqrt(n), 1/sqrt(n), ind, Llist);
        csvwrite(append(path,'Mcompare_',string(j),'.csv'),M(:,Llist));
        

    end
    
end
delete(gcp('nocreate'));
