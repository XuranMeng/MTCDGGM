function[betaoracle] = oracle(Z0, iU,beta)

[n,p]=size(Z0);
q=size(iU,2)-1;
dim=(p-1)*(q+1);
A=sparse(n*p,p*dim); %this is the large design matrix, which is sparse
b=[];
interM = [];

for j=1:(q+1)
    interM=[interM Z0.*iU(:,j)]; %these are the interaction terms between Z0 and iU
end


j=1;
beta_1=beta((1:(p-1)*(q+1))+(j-1)*(p-1)*(q+1),:);
Llist1=find(beta_1~=0);
X1=interM;
X1(:,((0:q)*p+j)) = [];
% A(((j-1)*n+1):(j*n),((j-1)*dim+1):(j*dim)) = Xj;
y1=Z0(:,j);
b1=y1-mean(y1);%-mean(y1)
X1oracle=X1(:,Llist1);




j=2;
beta_2=beta((1:(p-1)*(q+1))+(j-1)*(p-1)*(q+1),:);
Llist2=find(beta_2~=0);
X2=interM;
X2(:,((0:q)*p+j)) = [];
% A(((j-2)*n+2):(j*n),((j-2)*dim+2):(j*dim)) = Xj;
y2=Z0(:,j);
b2=y2-mean(y2);%-mean(y2)
X2oracle=X2(:,Llist2);



betaoracle=[inv(X1oracle' * X1oracle) * (X1oracle' * b1);inv(X2oracle' * X2oracle) * (X2oracle' * b2)];

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




