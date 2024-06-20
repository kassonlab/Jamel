import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
#Assumptions:1. Smoothing values should be tested on control and applied to all
#2. I can apply smoothing to first derivative
#3. Peaks that are less than half the max value are not important
#4 peaks should be far away from each other, in this case 100
#I need smoothed->smoother/rawderivatives and raw->rawderivative
#4 plots
peak_distance=100
def smooth(dataframe,window):
    return dataframe.rolling(window=window).mean()
def find_inflection_temperatures(data_column, dt, raw_smooth_window, derivative_smooth_window):
    smoothed_data = smooth(data_column, raw_smooth_window)
    raw_derivative_from_raw=data_column.diff()/dt
    raw_derivative_from_smooth=smoothed_data.diff()/dt
    smooth_first_derrivative_w_na=smooth(raw_derivative_from_smooth, derivative_smooth_window)
    peaks=find_peaks(smooth_first_derrivative_w_na,distance=peak_distance,height=smooth_first_derrivative_w_na.max()/2)[0]
    return smoothed_data,raw_derivative_from_smooth,smooth_first_derrivative_w_na,raw_derivative_from_raw,peaks

dsf_data=pd.read_csv(r"C:\Research\Total_DSF_30ug - Sheet1.csv")
# dsf_data.iloc[:,2:]=dsf_data.iloc[:,2:].sub(dsf_data['BALD'],axis=0)
# dsf_data.to_csv(r'C:\Research\dsf_data_subtracted.csv',index=False)
smoothed_data=dsf_data.copy(True)
raw_derivative_from_smoothed=dsf_data.copy(True)
smooth_derivative_from_smoothed=dsf_data.copy(True)
raw_derivative_from_raw=dsf_data.copy(True)
for x in range(2,11):
    test=dsf_data.iloc[:,x]
    column_name=dsf_data.columns[x]
    peaks=find_inflection_temperatures(test,0.2,40,40)
    smoothed_data.iloc[:, x]=peaks[0]
    raw_derivative_from_smoothed.iloc[:, x]=peaks[1]
    plt.plot(dsf_data.iloc[:, 1], raw_derivative_from_smoothed.iloc[:,x], label='Smoothed Data', color='red')
    smooth_derivative_from_smoothed.iloc[:, x]=peaks[2]
    plt.plot(dsf_data.iloc[:, 1], smooth_derivative_from_smoothed.iloc[:, x], label='Smoothed Data', color='blue')
    raw_derivative_from_raw.iloc[:, x]=peaks[3]
    plt.show()
# smoothed_data.to_csv(r'C:\Research\smoothed_data_subtracted.csv',index=False)
# smooth_derivative_from_smoothed.to_csv(r'C:\Research\smoothed_derivative_from_smooth_subtracted.csv',index=False)
# raw_derivative_from_smoothed.to_csv(r'C:\Research\raw_derivative_from_smoothed_subtracted.csv',index=False)
# raw_derivative_from_raw.to_csv(r'C:\Research\raw_derivative_from_raw_subtracted.csv',index=False)
# plt.figure()
# plt.plot(first_derrivative.iloc[:,1], first_derrivative['BALD_smooth'], label='Smoothed Data', color='red')
# plt.plot(dsf_data.iloc[:,1], dsf_data['BALD'], label='Original Data')
# plt.plot(smoothed_data.iloc[:,1], smoothed_data['BALD_smooth'], label='Smoothed Data', color='red')
# plt.legend()
# plt.show()