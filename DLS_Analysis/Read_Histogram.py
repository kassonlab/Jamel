import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def read_DLS_histo(csv_file):
    histogram=pd.read_csv(csv_file,skiprows=2,skipfooter=2,engine='python')
    # Getting rid of the second half DLS Histo than transitions into a different histogram
    word_index = histogram[histogram['Radius (nm)'] == 'Histogram'].index[0]
    histogram = histogram.iloc[:word_index]
    histogram=histogram.astype('float64')
    filtered_histo=histogram[histogram['Intensity']>0]
    return filtered_histo


def condense_intensity(diameter_array):
    dia_column=diameter_array[1:,0]
    inten_column=diameter_array[1:,1]
    diameter_range=max(dia_column)-min(dia_column)
    max_internsity=np.argmax(inten_column)
    prominent_diameter=dia_column[max_internsity]
    std_dev=np.std(dia_column)
    info_array=np.array([['Range','Max Intensity','High Inten. Diameter','Std Dev'],
                         [diameter_range,max_internsity,prominent_diameter,std_dev]],dtype=object)
    return info_array

