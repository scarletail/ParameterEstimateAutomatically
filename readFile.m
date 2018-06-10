function [] = readFile(filename)
%UNTITLED3 此处显示有关此函数的摘要
%   此处显示详细说明
load(filename);
sim('RC1_Manual','SrcWorkspace','current');
end