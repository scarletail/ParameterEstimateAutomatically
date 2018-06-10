import numpy as np
from matplotlib import pyplot as plt
import pandas as pd


filename = r'C:\Users\justs\PycharmProjects\ParametersAutoEstimater\-化成TEST-2018-04-20 08-41-56-CH3-1672\-化成TEST-2018-04-20 08-41-56-CH3-1672-detail.xls'
dataframe = pd.read_excel(filename)
vol = np.array(dataframe['电压(V)'])
cur = np.array(dataframe['电流(A)'])
tem = np.array(dataframe['温度(℃)'])
plt.figure()
plt.plot(vol)
plt.show()