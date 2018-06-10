import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


SOCinit = 9


def zeroState(t, a, b):
    return a*(1 - np.exp(-t/b))


def zeroInput(t, a, b):
    return a*np.exp(-t/b)


filename = r'_201706141616#001_1.xls'

data = pd.read_excel(filename, sheet_name='记录层').dropna(axis=0, how='any')
step = np.array(data['工步序号'])



# 寻找断点
breakpoint =np.array([0])
for i in range(1, len(step)):
    if step[i] != step[i - 1]:
        breakpoint = np.append(breakpoint, i)
breakpoint = np.append(breakpoint, len(step))
SOC = SOCinit

SOCList = np.array([])
UOCList = np.array([])
ROList1 = np.array([])
ROList2 = np.array([])
ROList3 = np.array([])
RList = np.array([])
CList = np.array([])
TaoList1 = np.array([])
TaoList2 = np.array([])
flag = 0   #开始标志
lastVol = 0
for i in range(len(breakpoint) - 1):
    subdata = data[breakpoint[i]: breakpoint[i + 1]]  #划分子块
    substatus = np.array(subdata['工步状态'])
    substep = np.array(subdata['工步序号'])[0]
    substep = int(substep)
    tim = np.array(subdata['相对时间(Min)'])
    vol = np.array(subdata['电压(V)'])
    cur = np.array(subdata['电流(A)'])

    if substatus[0] == '静置' and len(substatus) > 3000 :  #零输入响应，计算Tao
        flag = 1
        print('当前状态：静置；工步序号： ' + str(substep))
        print('SOC is ' + str(SOC))
        uoc = vol[-1]
        tim = tim[:100]
        vol = vol[:100]
        up = []
        for i in range(len(vol)):
            up.append( vol[i] - uoc )
        popt1, pcov1 = curve_fit(zeroInput, tim, up)
        b = popt1[1]
        print('b = ')
        print(b)
        TaoList2 = np.append(TaoList2, b)
    elif flag == 1:

        if substatus[0] == '恒流充电' and len(substatus) == 12:   #充电，无事发生
            print('当前状态：脉冲充电；工步序号：'+ str(substep))
            print('SOC is '+ str(SOC))
        elif substatus[0] == '恒流放电':  #SOC减0.1
            if len(substatus) <= 362 and len(substatus) >= 361:
                SOC -= 1
                print('当前状态：恒流放电（3Ah）；工步序号：' + str(substep))
                print('SOC is ' + str(SOC))
            elif len(substatus) == 12: #零状态响应，计算开路电压、欧姆内阻、R、C
                print('当前状态：脉冲放电；工步序号：' + str(substep))
                print('SOC is ' + str(SOC))
                r0 = (vol[1] - vol[0]) / cur[1]  #第一种R0
                ROList1 = np.append(ROList1, r0)
                uoc = vol[0]
                UOCList = np.append(UOCList, uoc)  #Uoc固定的计算方式
                f = 0
                c_sum = 0
                for i in range(1, len(cur)):
                    c_sum = c_sum + cur[i]
                    f = f + 1
                    c_avr = c_sum / f
                up = []
                for i in range(len(tim)):
                    temp = vol[i] - uoc - r0 * cur[i]
                    up.append(temp)
                popt1, pcov1 = curve_fit(zeroState, tim, up)
                a = popt1[0]
                b = popt1[1]
                R = a / c_avr
                C = b / R
                CList = np.append(CList, C)
                RList = np.append(RList, R)
                TaoList1 = np.append(TaoList1, b)
                SOCList = np.append(SOCList, round(SOC/10, 1))

TaoList2 = np.delete(TaoList2, -1)



df = pd.DataFrame({'SOC': SOCList,
                   'Uoc': UOCList,
                   'R0': ROList1,
                   'Tao from zero-state': TaoList1,
                   'Tao from zero-input': TaoList2,
                   'R': RList,
                   'C': CList})
columns = ['SOC','Uoc','R0','Tao from zero-state', 'Tao from zero-input', 'R','C']
df.to_csv("result.csv",index=False,sep=',',columns=columns)