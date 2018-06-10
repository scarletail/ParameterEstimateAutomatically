function [result] = dataProc(subdata)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
substatus = subdata.x_____3;
result = 'useless';
if strcmp(substatus{1}, '静置') && (length(substatus) > 3000)
    result = 'zeroinput';
end
if strcmp(substatus{1}, '恒流放电') && (length(substatus) <=362) && (length(substatus) >= 361)
    result = 'discharge';
end
if strcmp(substatus{1}, '恒流放电') && (length(substatus) == 12)
    result = 'zerostate';
end
end