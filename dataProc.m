function [result] = dataProc(subdata)
%UNTITLED3 �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
substatus = subdata.x_____3;
result = 'useless';
if strcmp(substatus{1}, '����') && (length(substatus) > 3000)
    result = 'zeroinput';
end
if strcmp(substatus{1}, '�����ŵ�') && (length(substatus) <=362) && (length(substatus) >= 361)
    result = 'discharge';
end
if strcmp(substatus{1}, '�����ŵ�') && (length(substatus) == 12)
    result = 'zerostate';
end
end