%%*************************************************************************
%% GMMReg:
%% Copyright (c) 2022 by
%% Jingfei Zhang and Yi Li
%%*************************************************************************
%% Input - Z0 centered responses of dimension n*p
%%       - iU n*q design matrix 
%% Note: alpha is fixed at 0.75 for faster computation
%%*************************************************************************

function[y_sepa,beta] = sepaGMMRegfixlbd(Z0, iU,lbde,lbdg)

[n,p]=size(Z0);
q=size(iU,2)-1;
dim=(p-1)*(q+1);
A=sparse(n*p,p*dim);
b=[];
interM = [];
for j=1:(q+1)
    interM=[interM Z0.*iU(:,j)];
end

%% group structure
int=ceil((1:dim)/(p-1));
G = [];
for i = 1:(q+1)
    G = [G;find(int==i)'];
end
ind = zeros(3,q+1);
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

%% nodewise regression
beta=zeros(p,(p-1)*(q+1));
y_sepa=zeros(p,n);
parfor j=1:p
    A=interM;
    A(:,((0:q)*p+j)) = [];
	y1=Z0(:,j);
	b=y1-mean(y1);
    lam1_max = norm(A'*b,Inf); 

    lambda1  = lam1_max*lbde;

    lambda2 = lam1_max*lbdg;
    Amap = @(x) mexMatvec(A,x,0);
    ATmap = @(y) mexMatvec(A,y,1);
    AATmap = @(x) Amap(ATmap(x));
    eigsopt = struct();
    eigsopt.issym = 1; 
    eigsopt.tol=1e-3;
    Lip = eigs(AATmap,length(b),1,'LA',eigsopt);

    c = [lambda1;lambda2];
    Amap = @(x) mexMatvec(A,x,0);
    ATmap = @(y) mexMatvec(A,y,1);
    AATmap = @(x) Amap(ATmap(x));
    opts = struct();
    opts.Lip = Lip;
    opts.stoptol = 1e-6;
    opts.printyes = 0;
    Ainput = struct();
    Ainput.A = A;
    Ainput.Amap = @(x) Amap(x);
    Ainput.ATmap = @(x) ATmap(x);
    [obj,y,z,x,info,runhist] = SGLasso_SSNAL(Ainput,b,size(A,2),c,G,ind,opts);
    beta(j,:)=x';
    y_sepa(j,:)=y';
end
end
