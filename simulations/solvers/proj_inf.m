function [y,rr] =  proj_inf(yinput,ld)
y = yinput;
normy = norm(yinput,Inf);
if normy > ld
    y = yinput/normy*ld;
end
rr = (yinput == y);
end

% y 是投影后的向量，其  L_\infty  范数不超过 ld。
% 	•	rr 是布尔向量，表示 yinput 中的元素是否在投影前已经满足条件。

% 例如，输入 [1, -5, 2, 4] 会被缩放为 [1, -3, 2, 3]，这样它的  L_\infty  范数正好等于 ld = 3。