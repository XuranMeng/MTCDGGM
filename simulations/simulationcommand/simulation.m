% cd /Users/emmazhang/Library/CloudStorage/OneDrive-EmoryUniversity/MtRegGMM/code/SparseGroupLasso-main;
clear; rng('default');
addpath(genpath(pwd));
parpool(5);

for h=1:50
filename = ['Z0', num2str(h),'.csv'];
Z0=csvread(filename,1,0);
filename = ['iU', num2str(h),'.csv'];
iU=csvread(filename,1,0);
[n,p]=size(Z0);
q=size(iU,2)-1;
dim=(p-1)*(q+1);
A=sparse(n*p,p*dim);
b=[];
interM = [];
for j=1:(q+1)
    interM=[interM Z0.*iU(:,j)];
end

for j=1:p
	Xj=interM;
    Xj(:,((0:q)*p+j)) = [];
    A(((j-1)*n+1):(j*n),((j-1)*dim+1):(j*dim)) = Xj;
	y1=Z0(:,j);
	b=[b;y1-mean(y1)];
end

%% Joint SGL
%% group structure
int=repmat(ceil((1:dim)/(p-1)),1,p);
G = [];
for i = 1:(q+1)
    G = [G;find(int==i)'];
end
ind = zeros(3,q+1);
grpsize = round(size(A,2)/(q+1)); 
for i = 1:(q+1)
    if i == 1
        ind(1,1) = 1; ind(2,1) = grpsize; ind(3,1) = 0;
    else
        ind(1,i) = ind(2,i-1) + 1;
        ind(2,i) = i*grpsize;
        ind(3,i) = 1;
    end
end

%% cross validation
lam1_max = norm(A'*b,Inf); alpha = 0.75; nl=50;
lambda1  = lam1_max*exp(linspace(log(1), log(0.2), nl));
lambda2 = lambda1*(1-alpha)/alpha;
Amap = @(x) mexMatvec(A,x,0);
ATmap = @(y) mexMatvec(A,y,1);
AATmap = @(x) Amap(ATmap(x));

eigsopt.issym = 1;
eigsopt.tol=1e-3;
Lip = eigs(AATmap,length(b),1,'LA',eigsopt);

cut = repelem(1:5,n/5); 
Cerror_SGL = zeros(5,nl); Cerror_L= zeros(5,nl); Cerror_G = zeros(5,nl);
for t=1:5
    foldidx1 = repmat(cut~=t,1,p);
    Atrain = A(foldidx1,:);
    btrain = b(foldidx1);
    foldidx2 = repmat(cut==t,1,p);
    Atest = A(foldidx2,:);
    btest = b(foldidx2);
    Amap = @(x) mexMatvec(Atrain,x,0);
    ATmap = @(y) mexMatvec(Atrain,y,1);
    AATmap = @(x) Amap(ATmap(x));
    opts.stoptol = 1e-4;
    opts.printyes = 0;
    opts.Lip = Lip;
    Ainput.A = Atrain;
    Ainput.Amap = @(x) Amap(x);
    Ainput.ATmap = @(x) ATmap(x);
    %% solver
    parfor i=1:nl
        c = [lambda1(i);lambda2(i)];
        [obj,y,z,x,info,runhist] = SGLasso_SSNAL(Ainput,btrain,size(Atrain,2),c,G,ind,opts);
        supp = find(x~=0);
        A_supp = Atrain(:,supp);
        x_supp = inv(A_supp'*A_supp)*A_supp'*btrain;
        x(supp) = x_supp;
        res = Atest*x - btest;
        Cerror_SGL(t,i) = Cerror_SGL(t,i)+sum(res.^2);  
        c1 = [lambda1(i);0];
        [obj,y,z,x,info,runhist] = SGLasso_SSNAL(Ainput,btrain,size(Atrain,2),c1,G,ind,opts);
        supp = find(x~=0);
        A_supp = Atrain(:,supp);
        x_supp = inv(A_supp'*A_supp)*A_supp'*btrain;
        x(supp) = x_supp;
        res = Atest*x - btest;
        Cerror_L(t,i) = Cerror_L(t,i)+sum(res.^2);  
        %c2 = [0;lambda2(i)*5]; %c2 = [0;lambda2(i)];
        %[obj,y,z,x,info,runhist] = SGLasso_SSNAL(Ainput,btrain,size(Atrain,2),c2,G,ind,opts);
        %supp = find(x~=0);
        %A_supp = Atrain(:,supp);
        %x_supp = inv(A_supp'*A_supp)*A_supp'*btrain;
        %x(supp) = x_supp;
        %res = Atest*x - btest;
        %Cerror_G(t,i) = Cerror_G(t,i)+sum(res.^2);  
    end
end
meanC = mean(Cerror_SGL,1); sdC = std(Cerror_SGL,1)/sqrt(5);
[minv,minid]=min(meanC);


%% estimation
c = [lambda1(minid);lambda2(minid)];
Amap = @(x) mexMatvec(A,x,0);
ATmap = @(y) mexMatvec(A,y,1);
AATmap = @(x) Amap(ATmap(x));
opts.Lip = Lip;
opts.stoptol = 1e-6;
Ainput.A = A;
Ainput.Amap = @(x) Amap(x);
Ainput.ATmap = @(x) ATmap(x);
[obj,y,z,x,info,runhist] = SGLasso_SSNAL(Ainput,b,size(A,2),c,G,ind,opts);
supp = find(x~=0);
A_supp = A(:,supp);
x_supp = inv(A_supp'*A_supp)*A_supp'*b;
x(supp) = x_supp;
filename = ['beta', num2str(h),'.csv'];
csvwrite(filename,x);


%% joint Lasso
meanC = mean(Cerror_L,1); sdC = std(Cerror_L,1)/sqrt(5);
[minv,minid]=min(meanC);
c = [lambda1(minid);0];
[obj,y,z,x,info,runhist] = SGLasso_SSNAL(Ainput,b,size(A,2),c,G,ind,opts);
supp = find(x~=0);
A_supp = A(:,supp);
x_supp = inv(A_supp'*A_supp)*A_supp'*b;
x(supp) = x_supp;
filename = ['bb', num2str(h),'.csv'];
csvwrite(filename,x);


%% joint Group Lasso
%meanC = mean(Cerror_G,1); sdC = std(Cerror_G,1)/sqrt(5);
%[minv,minid]=min(meanC);
%c = [0;lambda2(minid)*5]; %c = [0;lambda2(minid)];
%[obj,y,z,x,info,runhist] = SGLasso_SSNAL(Ainput,b,size(A,2),c,G,ind,opts);
%%supp = find(x~=0);
%A_supp = A(:,supp);
%x_supp = inv(A_supp'*A_supp)*A_supp'*b;
%x(supp) = x_supp;
%filename = ['bbb', num2str(h),'.csv'];
%csvwrite(filename,x);


%% Separate SGL
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

beta=[];
for j=1:p
    A=interM;
    A(:,((0:q)*p+j)) = [];
	y1=Z0(:,j);
	b=y1-mean(y1);
    lam1_max = norm(A'*b,Inf); 
    lambda1  = lam1_max*exp(linspace(log(1), log(0.2), nl));
    lambda2 = lambda1*(1-alpha)/alpha;
    cut = repelem(1:5,n/5); 
    Cerror = zeros(5,nl);
    Amap = @(x) mexMatvec(A,x,0);
    ATmap = @(y) mexMatvec(A,y,1);
    AATmap = @(x) Amap(ATmap(x));
    Lip = eigs(AATmap,length(b),1,'LA',eigsopt);

    for t=1:5
        foldidx1 = cut~=t;
        Atrain = A(foldidx1,:);
        btrain = b(foldidx1);
        foldidx2 = cut==t;
        Atest = A(foldidx2,:);
        btest = b(foldidx2);
        Amap = @(x) mexMatvec(Atrain,x,0);
        ATmap = @(y) mexMatvec(Atrain,y,1);
        AATmap = @(x) Amap(ATmap(x));
        opts.stoptol = 1e-4;
        opts.printyes = 0;
        opts.Lip = Lip;
        Ainput.A = Atrain;
        Ainput.Amap = @(x) Amap(x);
        Ainput.ATmap = @(x) ATmap(x);
        for i=1:nl
            c = [lambda1(i);lambda2(i)];
            [obj,y,z,x,info,runhist] = SGLasso_SSNAL(Ainput,btrain,size(Atrain,2),c,G,ind,opts);
            supp = find(x~=0);
            if length(supp)<size(Atrain,1)
                A_supp = Atrain(:,supp);
                x_supp = inv(A_supp'*A_supp)*A_supp'*btrain;
                x(supp) = x_supp;
            else
                x = x*100;
            end
            res = Atest*x - btest;
            Cerror(t,i) = Cerror(t,i)+sum(res.^2);  
        end
    end
    meanC = mean(Cerror,1); sdC = std(Cerror,1)/sqrt(5);
    [minv,minid]=min(meanC); 
   

    %% estimation
    c = [lambda1(minid);lambda2(minid)];
    Amap = @(x) mexMatvec(A,x,0);
    ATmap = @(y) mexMatvec(A,y,1);
    AATmap = @(x) Amap(ATmap(x));
    opts.Lip = Lip;
    opts.stoptol = 1e-6;
    Ainput.A = A;
    Ainput.Amap = @(x) Amap(x);
    Ainput.ATmap = @(x) ATmap(x);
    [obj,y,z,x,info,runhist] = SGLasso_SSNAL(Ainput,b,size(A,2),c,G,ind,opts);
    supp = find(x~=0);
    A_supp = A(:,supp);
    x_supp = inv(A_supp'*A_supp)*A_supp'*b;
    x(supp)= x_supp;
    beta=[beta,x'];
end
filename = ['b', num2str(h),'.csv'];
csvwrite(filename,beta');

end
