import os

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from os import listdir,path
temp_dict={'B':'25','C':'35','D':'55','E':'75'}
virus_dict={'6':'bat2006','7':'bat25','8':'bloom','2':'eid','5':'GCC','4':'sorex','3':'swine','1':'zhe'}
file_labels={}
for letter,temp in temp_dict.items():
    for number,virus in virus_dict.items():
        file_labels[letter+number]=f'{virus}_{temp}C'
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

def manage_DLS_files(file_label_dict:dict,DLS_directory:str):
    for file in [path.join(DLS_directory,file) for file in listdir(DLS_directory)]:
        for key,value in file_label_dict.items():
            if key in file:
                os.rename(file,file.replace(key,value))
                break



# radius_column='Radius (nm)'
# intensity_column='Intensity'
# particle_histo=pd.read_csv(r"C:\Research\DLS_data.csv")
# temp_dict={25: 'green', 35: 'blue', 55: 'pink', 75: 'red'}
# viruses=particle_histo['Virus'].unique()
# For creating a pseudo ridge plot where each temperature has its radius intensity stacked on top of each other as its own plot
# for virus in viruses:
#     virus_data=particle_histo[particle_histo['Virus']==virus]
#     fig, axes = plt.subplots(nrows=4, figsize=(10, 8))
#     temp_split_data = [(virus_data[virus_data['Temperature'] == x], x) for x in temp_dict.keys()]
#     for i, (data,temp) in enumerate(temp_split_data):  # Skip the first column 'Category'
#         axes[i].bar(data[radius_column],data[intensity_column])
#         axes[i].set_title(f'{temp}Celsius')
#         axes[i].set_xlabel(radius_column)
#         axes[i].set_ylabel(intensity_column)
#         axes[i].set_xlim(0, 150)
#     plt.title(virus)
#     plt.tight_layout()
#     plt.show()
    # plt.savefig(rf'C:\Research\{virus}_dls_ridgeplot.png')
# based of a 90% cutoff its 38.28 is good and everything else is bad
# we need bald DLS
# we need newer data with fresher samples

# good_radius_cutoff=38.28
# fraction_good_DLS:dict[str,list]={virus:[] for virus in viruses}
# For mapping Temperature with fraction of intensity under cutoff thats based on the bottom 90% of SARS data at 25Celsius
# for virus in viruses:
#     virus_data = particle_histo[particle_histo['Virus'] == virus]
#     temp_split_data = [(virus_data[virus_data['Temperature'] == x], x) for x in temp_dict.keys()]
#     plot_data={}
#     for i, (data, temp) in enumerate(temp_split_data):
#         good=data[data[radius_column]<=good_radius_cutoff]['Intensity'].sum()
#         bad=data[data[radius_column] > good_radius_cutoff]['Intensity'].sum()
#         fraction_good_DLS[virus].append(good/(good+bad) if (good+bad)>0 else 0)
# pd.DataFrame(fraction_good_DLS).to_csv(r'C:\Research\fraction_good_dls.csv',index=False)
#         plot_data[temp]=good/(good+bad) if (good+bad)>0 else 0
#     plt.plot(plot_data.keys(),plot_data.values())
# plt.legend(viruses)
# plt.xlabel('Temperature')
# plt.ylabel('Fraction Good')
# plt.show()
# Now each distribution is on top of each other per virus
# pltnum=1
# temp_split_data = [(particle_histo[particle_histo['Temperature'] == Cel], Cel,color) for Cel,color in temp_dict.items()]
# for temp_data,temp,color in temp_split_data:
#     plt.figure(pltnum)
#     pltnum+=1
#     for virus in viruses:
#         virus_data=temp_data[temp_data['Virus']==virus]
#         plt.plot(virus_data[radius_column],virus_data[intensity_column])
#         plt.xlim(0, 150)
#         plt.ylim(0,3)
#     plt.legend(viruses)
#     plt.title(temp)
#     plt.xlabel(radius_column)
#     plt.ylabel(intensity_column)
    # plt.show()
    # plt.savefig(rf'C:\Research\{temp}_overlayed_histo.png')
