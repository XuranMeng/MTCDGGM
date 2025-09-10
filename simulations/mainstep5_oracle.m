mex mexfun/mexbwsolve.c
mex mexfun/mexfwsolve.c
mex mexfun/mexfz.c
mex mexfun/mexMatvec.c
mex mexfun/mexProjL2.c
mex mexfun/mexProxL2.c
mex mexfun/mextriang.c
mex mexfun/mexMatvecmultiple.c
mex mexfun/mexmatrixmul.c
addpath('/Users/mengxuran/Documents/program/2024/MtGMMRegcode/solvers')
addpath('/Users/mengxuran/Documents/program/2024/MtGMMRegcode/tensor_toolbox')
addpath('/Users/mengxuran/Documents/program/2024/MtGMMRegcode/simulationcommand')
addpath('/Users/mengxuran/Documents/program/2024/MtGMMRegcode/else')
addpath('/Users/mengxuran/Documents/program/2024/MtGMMRegcode/SGL_functions')
savepath

% parfor i=1:10
%   i;
% end

nlist=[800 400 200 100]
% pairs = [20, 120; 120, 20];  % Define the pairs of (i, j)

Treptition=200;

for n=nlist
% for k = 1:size(pairs, 1)
p = 70;
q = 120;
    if n == 400
        startrep = 1;
    elseif n==200
        startrep=1
    else
        startrep = 1;
    end
% % Define parameters 

nu=0.3;
dim=(p-1)*(q+1);


path=append("datarepo/n=",string(n),"p=",string(p),"q=",string(q),"nu=",string(nu),'/');
beta=csvread(append(path,'beta.csv'));
% path=append("datarepo/n=",string(n),"p=",string(p),"q=",string(q),"nu=",string(nu),'/rep',string(1),'/');
% Z0=csvread(append(path,'Z0.csv'));
% iU=csvread(append(path,'iU.csv'));

% [y_mul,hatbeta_sepa]=sepaGMMReg(Z0,iU);
% hatbeta_sepa(abs(hatbeta_sepa)>0.10);
% a=reshape(hatbeta_sepa',[],1);
% a(abs(a)>0.00)

% rng(1118)
% generate(n,p,q,nu,Treptition)
oracle_bias=[];
for reptition=1:200
    
    disp([n p q reptition]);
    path=append("datarepo/n=",string(n),"p=",string(p),"q=",string(q),"nu=",string(nu),'/rep',string(reptition),'/');
    Z0=csvread(append(path,'Z0.csv'));
    iU=csvread(append(path,'iU.csv'));
    [hatbeta_oracle]=oracle(Z0,iU,beta);
    oracle_bias=[oracle_bias;hatbeta_oracle'+0.3];
    % find(hatbeta_mul~=0)'
    % hatbeta_mul(abs(hatbeta_mul)>0.1)
    % csvwrite(append(path,'hatbeta_muloracle.csv'),hatbeta_mul);
    % csvwrite(append(path,'y_mul.csv'),y_mul);
    % [y_sepa,hatbeta_sepa]=sepaGMMReg(Z0,iU);
    % hatbeta_sepa(hatbeta_sepa~=0)'
    % a=reshape(hatbeta_sepa',[],1);
    % csvwrite(append(path,'hatbeta_sepa.csv'),a);
    % csvwrite(append(path,'y_sepa.csv'),y_sepa);
    % 
end
path=append("datarepo/n=",string(n),"p=",string(p),"q=",string(q),"nu=",string(nu),'/');
csvwrite(append(path,'oracle_bias.csv'),oracle_bias);
% end
end