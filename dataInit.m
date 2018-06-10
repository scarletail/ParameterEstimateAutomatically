function [Capacity, CurrentA, Qe_init, VoltageV] = dataInit(discFile)
Capacity = 31;
Qe_init = 0;
discdata = readtable(discFile);
discdata(1,:) = [];
CurrentA = discdata.x___1;
VoltageV = discdata.x___2;
end