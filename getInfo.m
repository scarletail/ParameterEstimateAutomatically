function [t, I_in, V_out] = getInfo(newdata)
%UNTITLED9 此处显示有关此函数的摘要
%   此处显示详细说明
I_in = newdata.x___A_;
t = 1:height(newdata);
t = t';
V_out = newdata.x___V_;
end
