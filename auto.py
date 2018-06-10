import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
def func1(t, a, b):
    return a*(1 - np.exp(-t/b))
# def func2(t, a, b, c, d):
#     return a*(1 - np.exp(-t/b)) + c*(1 - np.exp(-t/d))
filename = r'data.xlsx'
data = pd.read_excel(filename)
tim = data['相对时间(Min)']
vol = data[ '电压(V)']
cur = data[ '电流(A)']
uoc = vol[0]
r0 = (vol[1] - vol[0])/cur[1]
up=[]
c_sum = 0
flag = 0
for i in range(1, len(cur)):
    c_sum = c_sum + cur[i]
    flag = flag + 1
    c_avr = c_sum/flag
for i in range(len(tim)):
    temp = vol[i] - uoc - r0 * cur[i]
    up.append(temp)
popt1, pcov1=curve_fit(func1, tim, up)
a = popt1[0]
b = popt1[1]
R = a/c_avr
C = b/R
print('R = '+ str(R))
print('C = '+ str(C))
# popt2, pcov2=curve_fit(func2, tim, up)
# a = popt2[0]
# b = popt2[1]
# c = popt2[2]
# d = popt2[3]
# R1 = a/c_avr
# C1 = b/R1
# R2 = c/c_avr
# C2 = d/R2
# print('R1 = '+ str(R1))
# print('C1 = '+ str(C1))
# print('R2 = '+ str(R2))
# print('C2 = '+ str(C2))