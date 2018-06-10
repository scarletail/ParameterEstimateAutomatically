import pandas as pd
import numpy as np
from scipy.optimize import curve_fit

def zeroInput(t, a, b):
    return a*np.exp(-t/b)

filename =r'analysis.xlsx'
df = pd.read_excel(filename)
up = np.array(df['up'])
t = np.array(df['t'])
popt1, pcov1 = curve_fit(zeroInput, t, up)
a = popt1[0]
b = popt1[1]

print(a)
print(b)