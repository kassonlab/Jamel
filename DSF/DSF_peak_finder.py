import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
def smooth(dataframe,window):
    return dataframe.rolling(window=window).mean()
def find_inflection_temperatures(data_column,dt,raw_smooth,derivative_smooth):
    smoothed_data = smooth(data_column, raw_smooth)
    smoothed_data = smoothed_data.dropna()
    first_derrivative=smoothed_data.diff()/dt
    smooth_first_derrivative=smooth(first_derrivative,derivative_smooth)
    smooth_first_derrivative=smooth_first_derrivative.dropna()
    peaks=find_peaks(smooth_first_derrivative,distance=100,height=smooth_first_derrivative.max()/2)[0]
    return peaks

dsf_data=pd.read_csv(r"C:\Research\Total_DSF_30ug - Sheet1.csv")

# plt.plot(first_derrivative.iloc[:,1], first_derrivative['BALD_smooth'], label='Original Data')
for x in range(2,11):
    test=dsf_data.iloc[:,x]
    print(dsf_data.columns[x])
    peaks=find_inflection_temperatures(test,0.2,40,50)
    print(dsf_data.iloc[:,1].iloc[peaks])
# plt.figure()
# plt.plot(first_derrivative.iloc[:,1], first_derrivative['BALD_smooth'], label='Smoothed Data', color='red')
# plt.plot(dsf_data.iloc[:,1], dsf_data['BALD'], label='Original Data')
# plt.plot(smoothed_data.iloc[:,1], smoothed_data['BALD_smooth'], label='Smoothed Data', color='red')
# plt.legend()
# plt.show()