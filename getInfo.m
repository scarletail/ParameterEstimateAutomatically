function [t, I_in, V_out] = getInfo(newdata)
%UNTITLED9 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
I_in = newdata.x___A_;
t = 1:height(newdata);
t = t';
V_out = newdata.x___V_;
end
