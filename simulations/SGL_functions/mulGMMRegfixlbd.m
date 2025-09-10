function[y,beta] = mulGMMRegfixlbd(Z0, iU,lbde,lbdg)

[n,p]=size(Z0);
q=size(iU,2)-1;
dim=(p-1)*(q+1);
A=sparse(n*p,p*dim); %this is the large design matrix, which is sparse
b=[];
interM = [];
for j=1:(q+1)
    interM=[interM Z0.*iU(:,j)]; %these are the interaction terms between Z0 and iU
end

%% populates the design matrix with nonzero blocks

for j=1:p 
	Xj=interM;
    Xj(:,((0:q)*p+j)) = [];
    A(((j-1)*n+1):(j*n),((j-1)*dim+1):(j*dim)) = Xj;
	y1=Z0(:,j);
	b=[b;y1-mean(y1)];%-mean(y1)
end
% size(A)
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

%% 5-fold cross validation
lam1_max = norm(A'*b,Inf); 
lambda1  = lam1_max*lbde; %try 10 lambda1
lambda2 = lam1_max*lbdg; %we set alpha and do not tune lambda2


%% estimation
c = [lambda1;lambda2];
% c = [lam1_max*0.9;lam1_max*0.9/3];
Amap = @(x) mexMatvec(A,x,0);
ATmap = @(y) mexMatvec(A,y,1);
AATmap = @(x) Amap(ATmap(x));
eigsopt.issym = 1;
eigsopt.tol=1e-3;
Lip = eigs(AATmap,length(b),1,'LA',eigsopt);

opts.Lip = Lip;
opts.stoptol = 1e-4;
Ainput.A = A;
Ainput.Amap = @(x) Amap(x);
Ainput.ATmap = @(x) ATmap(x);
[obj,y,z,x,info,runhist] = SGLasso_SSNAL(Ainput,b,size(A,2),c,G,ind,opts);
% supp = find(x~=0);
% A_supp = A(:,supp);
% x_supp = inv(A_supp'*A_supp)*A_supp'*b; %we do a refitting here after selection to reduce bias from lasso
% x(supp) = x_supp;
beta=x;
end
% a=x(supp);
% length(supp)
% a(abs(a)>0.1)

% j=1
% Xj=interM;
% Xj(:,((0:q)*p+j)) = [];
% A(((j-1)*n+1):(j*n),((j-1)*dim+1):(j*dim)) = Xj;
% A(1:n,1:(p-1)*(q+1))


% csvwrite('beta.csv',x);




