function [r1,tao1] = zerostate1(t,up)
%UNTITLED6 此处显示有关此函数的摘要
%   此处显示详细说明
zs1 = @(a, t)-a(1)*(1-exp(-t/a(2)));  
%function handle zs1 is used for the first order zero-state model
zs2 = @(a, t)-a(1)*(1-exp(-t/a(2)))-a(3)*(1-exp(-t/a(4)));
%function handle zs2 is used for the second order zero-state model
init0 = [0.0001, 0.01];
init1 = [0.1, 0.1, 0.1, 0.1];
result1 = nlinfit(t, up, zs1, init0);
r1 = result1(1);
tao1 = result1(2);
end
