function [k,xnew] = cardcal(x,r)
n = length(x);
normx1 = norm(x,1);
[absx,idx] = sort(abs(x),'descend');
for i = 1:n
    if sum(absx(1:i)) >= r*normx1
       k = i;
       break;
    end
end
if nargout > 1
   xnew = zeros(n,1);
   idxnew = idx(1:k);
   xnew(idxnew) = x(idxnew);
end

% 输入：

% 	•	x: 一个向量。
% 	•	r: 一个比例，用于确定稀疏度。

% 输出：

% 	•	k: 向量的稀疏度（即有多少个元素的绝对值和达到给定比例 r）。
% 	•	xnew: 稀疏化后的向量，保留前 k 个最大元素，其他位置为零（如果需要这个输出）。
