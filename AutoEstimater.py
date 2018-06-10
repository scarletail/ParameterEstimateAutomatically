import numpy as np
import pandas as pd
from scipy.optimize import curve_fit


SOCinit = 0.9
fInit =8
def zeroInput(t, a, b):
    return a*(1 - np.exp(-t/b))


def zeroStatus(t, a, b):
    return a*np.exp(-t/b)


fname = r'datasource.xlsx'
data = pd.read_excel(fname).dropna(axis=0,how='any')
step = np.array(data['工步序号'])


# 寻找断点
breakpoint =np.array([0])
for i in range(1, len(step)):
    if step[i] != step[i - 1]:
        breakpoint = np.append(breakpoint, i)
breakpoint = np.append(breakpoint, len(step))
SOC = SOCinit
f = fInit
SOCList = np.linspace(0.1,0.9,9)
UOCList = np.array([0.0] * 9)
ROList = np.array([0.0]*9)
RList = np.array([])
CList = np.array([])
for i in range( len(breakpoint) - 1 ):
    subdata = data[breakpoint[i]: breakpoint[i + 1]]  #划分子块
    substatus = np.array(subdata['工步状态'])
    substep = np.array(subdata['工步序号'])[0]
    substep = int(substep)

    if substatus[0] == '恒流充电':   #充电，无事发生
        print('当前状态：恒流充电；工步序号：'+ str(substep))
        print('SOC is '+ str(SOC))
    elif substatus[0] == '恒流放电':  #SOC减0.1
        if len(substatus) > 12:
            SOC = round(SOC - 0.1, 1)
            f = f - 1
            print('当前状态：恒流放电（3Ah）；工步序号：' + str(substep))
            print('SOC is ' + str(SOC))
        else: #零状态响应，计算开路电压、欧姆内阻、R、C
            print('当前状态：恒流放电；工步序号：' + str(substep))
            print('SOC is ' + str(SOC))
            tim = np.array(subdata['相对时间(Min)'])
            vol = np.array(subdata['电压(V)'])
            cur = np.array(subdata['电流(A)'])
            uoc = vol[0]
            UOCList[f] = uoc
            r0 = (vol[1] - vol[0]) / cur[1]
            ROList[f] = r0
            up = []
            c_sum = 0
            flag = 0
            for i in range(1, len(cur)):
                c_sum = c_sum + cur[i]
                flag = flag + 1
                c_avr = c_sum / flag
            for i in range(len(tim)):
                temp = vol[i] - uoc - r0 * cur[i]
                up.append(temp)
            popt1, pcov1 = curve_fit(zeroInput, tim, up)
            a = popt1[0]
            b = popt1[1]
            R = a / c_avr
            C = b / R
            RList = np.append(RList, R)
            CList = np.append(CList, C)
            print('R = ' + str(R))
            print('C = ' + str(C))

    else:  #零输入响应，忽略欧姆内阻计算R、C
        print('当前状态：静置；工步序号： ' + str(substep))
        print('SOC is ' + str(SOC))
df = pd.DataFrame({'SOC': SOCList[::-1], 'Uoc':UOCList[::-1], 'R0':ROList[::-1], 'R':RList, 'C':CList})
columns = ['SOC','Uoc','R0','R','C']
df.to_csv("result.csv",index=False,sep=',',columns=columns)