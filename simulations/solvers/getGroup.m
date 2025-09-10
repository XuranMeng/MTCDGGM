function [newA,G,ind,reorder] = getGroup(A,grpNUM)
% A,b      - data
% grpNUM   - number of groups

% G,ind    - input for Def_P.m(weight wi = sqrt(|G_i|))
[p] = size(A,2); %p:number of variables
if grpNUM > p, error('grpNUM > #variables'); end
ind = zeros(3,grpNUM);

% •	ind 是一个 3 × grpNUM 的矩阵，用来存储每个分组的起始索引、终止索引和该分组的大小。
% •	G 是从 1 到 p 的索引向量。

G = [1:p];

%% reordering p variables
reorder = randperm(p);
newA = A(:,reorder);

%% partition [1:p] into grpNUM groups

% probs = rand(nG,1); probs = probs / sum(probs); 
% probs = ones(nG,1)/nG;

probs = 1 + 0.3 * sign(randn(grpNUM,1)) .* rand(grpNUM,1); 
probs = probs / sum(probs); 
probs = cumsum(probs);

for i = 1:grpNUM
    if i == 1
        tmp = round(probs(1)*p); tmp = max(tmp,1);
        ind(1,1) = 1; ind(2,1) = tmp; ind(3,1) = sqrt(tmp);
    else
        ind(1,i) = ind(2,i-1) + 1;
        ind(2,i) = max(round(probs(i)*p),ind(1,i));
        ind(3,i) = sqrt(ind(2,i)-ind(1,i));
    end
end
% •	第一个分组的起始索引为 1，终止索引通过 probs 计算并四舍五入到最近的整数，大小的平方根存储在 ind(3,1)。
% •	后续的分组起始索引为上一个分组的终止索引的下一个位置，终止索引通过 probs 计算，分组大小的平方根存储在 ind 的第 3 行中。

end


% 这段代码定义了一个函数 getGroup，用于将数据矩阵的列（变量）重新分组，并返回分组信息和经过重新排列的矩阵。

% 输入：

% 	•	A: 数据矩阵，大小为 n x p，其中 n 是样本数，p 是变量数。
% 	•	grpNUM: 分组的数量。

% 输出：

% 	•	newA: 重新排列列后的数据矩阵。
% 	•	G: 变量的索引向量。
% 	•	ind: 包含每个分组的起始、终止索引和分组长度的矩阵。
% 	•	reorder: 重新排列后的变量顺序。

% 简要说明：

% 	1.	随机重新排列矩阵 A 的列（即变量），并将 p 个变量划分为 grpNUM 个组。
% 	2.	每组的大小根据随机生成的概率进行划分，并存储在 ind 中，ind 包括每组的起始和终止位置，以及对应分组大小的平方根。


% n = 12; % 样本数量
% p = 9;  % 变量数量
% grpNUM = 3; % 分组数量

% % 随机生成 12x9 的矩阵 A
% A = randn(n, p);

% % 调用 getGroup 函数
% [newA, G, ind, reorder] = getGroup(A, grpNUM);

% % 展示 ind 矩阵
% disp('ind 矩阵:');
% disp(ind);

% % 展示 G 向量
% disp('G 向量:');
% disp(G);

% % 展示 reorder 向量
% disp('reorder 向量:');
% disp(reorder);

% % 展示新的矩阵 newA
% disp('重新排列后的矩阵 newA:');
% disp(newA);