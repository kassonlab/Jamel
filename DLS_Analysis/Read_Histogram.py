import os
from collections.abc import Iterable
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from os import listdir, path

from matplotlib.ticker import MaxNLocator

temp_dict = {'B': '25', 'C': '35', 'D': '55', 'E': '75'}
virus_dict = {'6': 'bat2006', '7': 'bat25', '8': 'bloom', '2': 'eid', '5': 'GCC', '4': 'sorex', '3': 'swine',
              '1': 'zhe'}
file_labels = {}
for letter, temp in temp_dict.items():
    for number, virus in virus_dict.items():
        file_labels[letter + number] = f'{virus}_{temp}C'


def read_DLS_histo(csv_file):
    histogram = pd.read_csv(csv_file, skiprows=2, skipfooter=2, engine='python')
    # Getting rid of the second half DLS Histo than transitions into a different histogram
    word_index = histogram[histogram['Radius (nm)'] == 'Histogram'].index[0]
    histogram = histogram.iloc[:word_index]
    histogram = histogram.astype('float64')
    filtered_histo = histogram[histogram['Intensity'] > 0]
    return filtered_histo


def condense_intensity(diameter_array):
    dia_column = diameter_array[1:, 0]
    inten_column = diameter_array[1:, 1]
    diameter_range = max(dia_column) - min(dia_column)
    max_internsity = np.argmax(inten_column)
    prominent_diameter = dia_column[max_internsity]
    std_dev = np.std(dia_column)
    info_array = np.array([['Range', 'Max Intensity', 'High Inten. Diameter', 'Std Dev'],
                           [diameter_range, max_internsity, prominent_diameter, std_dev]], dtype=object)
    return info_array


def manage_DLS_files(file_label_dict: dict, DLS_directory: str):
    for file in [path.join(DLS_directory, file) for file in listdir(DLS_directory)]:
        for key, value in file_label_dict.items():
            if key in file:
                os.rename(file, file.replace(key, value))
                break


temp_dict = {25: 'green', 35: 'blue', 55: 'pink', 75: 'red'}


# (38.28 is good and everything else is bad
def stacked_histo_plots(dataframe: pd.DataFrame, subplot_variable_name: str, x_column='Radius (nm)',
                        y_column='Intensity', save_fig=''):
    # For creating a pseudo ridge plot where each temperature has its radius intensity stacked on top of each other as its own plot
    subplot_variables = dataframe[subplot_variable_name].unique()
    fig, axes = plt.subplots(nrows=len(subplot_variables), figsize=(10, 8))
    split_data = [(dataframe[dataframe[subplot_variable_name] == x], x) for x in subplot_variables]
    for i, (series, variable) in enumerate(split_data):  # Skip the first column 'Category'
        axes[i].bar(series[x_column], series[y_column])
        axes[i].set_title(f'{variable}')
        axes[i].set_xlabel(x_column)
        axes[i].set_ylabel(y_column)
        # axes[i].set_xlim(0, 150)
    plt.title(virus)
    plt.tight_layout()
    plt.show()
    plt.savefig(save_fig)


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
def DLS_summary_to_pandas(csv_file: str) -> pd.DataFrame:
    dataframe = pd.read_csv(csv_file, skiprows=1)
    return dataframe.iloc[1:]


def overlapping_line_plot(dataframes: list[pd.DataFrame] | tuple[pd.DataFrame], legend: list[str], x_column: str,
                          y_column: str, save_fig: str = ''):

    for plt_num,dataframe in enumerate(dataframes):
        # plt.figure(plt_num)
        plt.plot(dataframe[x_column], dataframe[y_column])
        plt.legend(legend)
        plt.xlabel(x_column)
        plt.ylabel(y_column)
    plt.gca().xaxis.set_major_locator(MaxNLocator(nbins=5))
    plt.savefig(save_fig) if save_fig else plt.show()
    plt.clf()


def find_t_half_index(series: pd.Series) -> int:
    half_value = (series.iloc[0]+series.iloc[-1])/2
    indexes_below_half = series.loc[series <= half_value].index
    return indexes_below_half[0]
# raw_summary files
if __name__ == '__main__':
    DLS_direc = r"C:\Research\DLS\ACF"
    # manage_DLS_files(file_labels,DLS_direc)
    ACF_files = [path.join(DLS_direc, file) for file in os.listdir(DLS_direc) if 'ACF' in file]
    for virus in virus_dict.values():
        virus_file = tuple(file for file in ACF_files if virus in file)
        virus_data = [DLS_summary_to_pandas(csv) for csv in virus_file]
        overlapping_line_plot(virus_data,[file.split("\\")[-1] for file in virus_file],'Time (s)',' Intensity')


# if __name__ == '__main__':
#     DLS_direc = r"C:\Research\DLS"
#     ACF_files = [path.join(DLS_direc, file) for file in os.listdir(DLS_direc) if 'total' in file]
#     for virus in virus_dict.values():
#         virus_file = tuple(file for file in ACF_files if virus in file)
#         virus_data = pd.read_csv(virus_file[0])
#         t_half_series=virus_data.apply(lambda x: virus_data.loc[find_t_half_index(x),'Time (s)'])[1:]
#         t_half_series.index=tuple(int(index.split('_')[-1][:2]) for index in t_half_series.index)
#         # t_half_series = t_half_series[~t_half_series.index.duplicated(keep='first')].sort_index()
#         t_half_series.plot()
#     plt.xlabel('Temperature')
#     plt.ylabel("Half-life (s)")
#     plt.legend(virus_dict.values())
#     plt.title('Half-lifes')
#     # plt.show()
#     plt.savefig(fr"C:\Research\DLS\All_ACF_half_life.png")

