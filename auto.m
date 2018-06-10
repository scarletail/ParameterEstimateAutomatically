filename = '_201706141616#001_1.xls';
discFile = '5C放电数据(1).xlsx';
[status,sheets] = xlsfinfo(filename);
%choose sheet
i = 5;
data = rmmissing(readtable(filename, 'sheet', sheets{i}));  %drop nan
%divide blocks
div = 0;
SOC = 1;
SOC_LUT = [];
R0 = [];
R1 = [];
C1 = []; 
Em = [];
newdata = cell2table({});
stepnum = data.x_____1;
for j = 1:(height(data) - 1)
    if stepnum(j + 1) ~= stepnum(j)
        div = [div, j];
    end
end
div = [div, height(data)];
flag = 0;
cond = 1;
for m = 1:(length(div)-1)
    subdata = data( (div(m)+1) : div(m+1) , : );
    %process subdata
    result = dataProc(subdata);
    if cond
        if strcmp(result, 'discharge')
            flag = 1;
            SOC = SOC - 0.1;
            SOC_LUT = [SOC_LUT, SOC];
            
        end
        if flag
            if strcmp(result, 'zeroinput')
                %change of flag mean to begin the whole process
                %begin zeroinput
                subdata = subdata(1:100, :);
                if abs(SOC) < 0.0001
                    cond = 0;
                end
            end
            
            if strcmp(result, 'zerostate')
                %begin zerostate
                [t, vol, cur] = getdata(subdata);
                r0 = ones(size(t)).*abs((vol(2)-vol(1))/cur(2));
                R0 = [R0; r0(1)];
                Em = [Em; vol(1)];
                uoc = ones(size(t)).*vol(1);
                up = vol - uoc - r0.*cur;
                [r1c, tao1] = zerostate1(t, up);
                nzc = nonzeros(cur);
                c_avr = sum(nzc(:))/length(nzc);
                r1 = abs(r1c/c_avr);
                c1 = tao1/r1;
                R1 = [R1; r1];
                C1 = [C1; c1];
                %checked, all right!
            end
            newdata = vertcat(newdata, subdata);
        end
    end
end
[t, I_in, V_out] = getInfo(newdata);
[Capacity, CurrentA, Qe_init, VoltageV] = dataInit(discFile);
SOC_LUT = SOC_LUT(:,1:length(R1));
sim('RC1_Manual');
%plot(t, V_out);
%hold on;
%plot(ScopeData1.time,ScopeData1.signals.values(:,1));