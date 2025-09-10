

addpath(genpath(pwd));

 % sample size
n=178;
p = 73; % dimension of response
q = 120; % dimension of covariates
nu=0.3;
dim=(p-1)*(q+1);


grpsize = p-1; 
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



hatbeta_mul=csvread('betaandM/hatbeta_mul.csv');
Llist=find(hatbeta_mul~=0);





Z0=csvread('Z0.csv',1,0);
iU=csvread('iU.csv',1,0);

interM = [];
for j=1:(q+1)
interM=[interM Z0.*iU(:,j)]; %these are the interaction terms between Z0 and iU
end


for j=1:p
    disp(j);
    betamul_j = hatbeta_mul((1:(p-1)*(q+1)) + (j-1)*(p-1)*(q+1), :);
    Llistmul = find(betamul_j ~= 0);
    Wj=interM;
    Wj(:,((0:q)*p+j)) = [];

    % Sigma_hat=mexmatrixmul(Wj',Wj)/n;
    [n,dim]=size(Wj);
    Xi=mexmatrixmul(Wj,Wj')/n;
    [U,D]=eig(Xi);
    if ~isempty(Llistmul)
        csvwrite(append('betaandM/', 'Mmulindex_', string(j), '.csv'), Llistmul);
        Mmul = OptimizeMnew(U, D, Wj, 2/sqrt(n), 1/sqrt(n), ind, Llistmul);
        csvwrite(append('betaandM/', 'Mmul_', string(j), '.csv'), Mmul(:,Llistmul));
    end

 
end

 
