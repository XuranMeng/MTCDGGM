
addpath(pwd)
n=200
p = 20; % dimension of response
q = 20; % dimension of covariates
nu=0.3;
dim=(p-1)*(q+1);
   
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

   
%% true coefficient represented as a p*p*(q+1) array
tB = tensor(zeros(1,p*p*(q+1)), [p p q+1]); 

rng(42);
listi=2:3;
listj=[1];
for i=listi
	for j=listj
		tB(j,j+1,i) = nu;
        tB(j+1,j,i) = tB(j,j+1,i);
    end
end

N = randi([0 1],n,q);
U = (N-0.5)/0.5;  %对列做归一化
iU = [ones(n,1), U];
X = zeros(n, p);
for i=1:n
    omega = ttv(tB,iU(i,:)', 3);
    omega = double(omega);
    omega = omega + 1.0*diag(ones(1,p));
    sigma = inv(omega);
    X(i,:) = mvnrnd(zeros(1,p),sigma,1);
end

Z0=X;


interM = [];
for j=1:(q+1)
    interM=[interM Z0.*iU(:,j)]; %these are the interaction terms between Z0 and iU
end



j=1;
% Llist=[1;2;Llist];
Wj=interM;
Wj(:,((0:q)*p+j)) = [];
Xi=mexmatrixmul(Wj,Wj')/n;
Sigma_hat=mexmatrixmul(Wj',Wj)/n;


