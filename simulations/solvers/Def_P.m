%% Define P map. 
%% Zhang Yangjing, 05 May 2017
%
%
% Input - prescribed group G_i, i=1,...,g
%       - G\in R^s, a column vector containing the indices of all the 
%         overlapping groups. G =[ G_1; G_2; ...; G_g]
%       - ind \in R^{3*g}, indices, G(ind(1,i):ind(2,i))denotes the indices
%         for the i-th group,
%         ind(3,i)= wi denotes the weight for the i-th group, 
%         e.g., we can take (wi)^2 = |G_i|.
%
% Output - P
%
% P.matrix   - matrix representation of P  
%
% P.grpNUM   - number of groups
% P.grpLEN   - a column vector [|G_1|, |G_2|, ... ,|G_grpNUM|]'
% P.ntotal   - |G_1|+|G_2|+...+|G_grpNUM|
%
%
% P.times    - P.times(x) computes P(x)%active
% P.trans    - P.trans(y) computes P'(y)%active
%
% P.Pi       - P.Pi(i) matrix representation of Pi, Pi(x)=x_{G_i}
%
%
% P.ProjL2   - P.ProjL2(z,c2) computes projection onto B_2 = B1*B2*...
%              Bi={zi \in R^|G_i|: norm(zi,2)<=c2*wi},i=1,2,...,grpNUM %active
%
% P.Lasso_fx - P.Lasso_fx(x) computes sum_i {wi*norm(P_i(x),2)}
% P.Lasso_fz - P.Lasso_fz(z) computes sum_i {wi*norm(zi,2)}

function P = Def_P(n,G,ind)

grpNUM = size(ind,2);
P.grpNUM = grpNUM;
grpLEN = (ind(2,:) - ind(1,:))' + 1;
P.grpLEN = grpLEN;

ntotal = sum(grpLEN);
P.ntotal = ntotal;
P.ind = ind;
P.G = G;

I = [1:ntotal];
J = G';
V = ones(1,ntotal);


Pma = sparse(I,J,V);
P.matrix = Pma;
 
P.Pi = @(i) P_i(i,G,ind,n,grpLEN);

P.times = @(x) Pma*x;
P.trans = @(y) Pma'*y;

P.ProjL2 = @(z,c1) mexProjL2(z,c1,ind,grpNUM);
P.ProxL2 = @(z,c1) mexProxL2(z,c1,ind,grpNUM);

P.Lasso_fx = @(x)mexfz(Pma*x,ind,grpNUM);
P.Lasso_fz = @(z)mexfz(z,ind,grpNUM);


end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Pi] = P_i(i,G,ind,n,grpLEN)
tmp = grpLEN(i);
I = [1:tmp];
J = G(ind(1,i):ind(2,i));
V = ones(1,tmp);
Pi = sparse(I,J,V,tmp,n);
end



% 输入：

% 	•	n: 矩阵的总列数。
% 	•	G: 一个用于分组的向量。
% 	•	ind: 用于定义每个分组的起止索引的矩阵。

% 输出：

% 	•	P: 一个结构体，包含以下字段：
% 	•	grpNUM: 分组的数量。
% 	•	grpLEN: 每个分组的长度。
% 	•	ntotal: 所有分组元素的总数。
% 	•	ind: 输入的分组索引。
% 	•	G: 输入的分组向量。
% 	•	matrix: 构造的稀疏矩阵。
% 	•	Pi: 一个函数句柄，用于提取第 i 个分组的稀疏矩阵。
% 	•	times: 一个函数句柄，表示稀疏矩阵与向量的乘法。
% 	•	trans: 一个函数句柄，表示稀疏矩阵的转置与向量的乘法。
% 	•	ProjL2 和 ProxL2: 函数句柄，表示一些 L2 投影和近似的操作。
% 	•	Lasso_fx 和 Lasso_fz: 用于处理 Lasso 问题的函数。


% G = [1, 2, 3, 2, 3, 4, 5, 6, 7];
% ind = [1, 4, 8; 3, 7, 9; sqrt(3), sqrt(4), sqrt(2)];
% n = 10;

% (1,1)        1
% (2,2)        1
% (3,3)        1

% (4,2)        1
% (5,3)        1
% (6,4)        1
% (7,5)        1

% (8,6)        1
% (9,7)        1
