import os
import pandas as pd
import statsmodels.api as sm
from statsmodels.formula.api import ols
import statsmodels.stats.multicomp as mc
from matplotlib import pyplot as plt
import numpy as np
import random
import DLS_Analysis.Read_Histogram
from scipy.stats import levene


def cum_dis_func(series_data: pd.Series,percentile_cutoff:float=0.9):
    x = series_data.astype(float).sort_values()
    y = np.arange(len(x)) / float(len(x))
    plt.plot(x, y)
    plt.show()
    for index,num in enumerate(y):
        if percentile_cutoff<num:
            return x[index]


def flatten_microscopy_files(grand_parent_directory, target_directory):
    """Dicrectories must be in specific structure 'grandparent/virus/temperature/fov/'"""
    for root, dirs, files in os.walk(grand_parent_directory):
        if not dirs:
            file = next((file for file in files if file.endswith('tif')))
            old_path = os.path.join(root, file)
            virus, temperature, fov = root.split('\\')[-3:]
            temperature = temperature[0:2]
            virus = virus.split('_')[-1]
            new_path = os.path.join(target_directory, f'{virus}_{temperature}_{fov}.ome.tif')
            os.system(f'copy {old_path} {new_path}')


# im trying a two way anova
def filter_temp_from_particle_table(particle_data: pd.DataFrame, temperature: int | str):
    particle_data=particle_data[particle_data['Temperature']==temperature]
    return particle_data.groupby('Virus')



# Performing two-way ANOVA
def two_way_anova(dataframe: pd.DataFrame, variable_one, variable_2, result_column):
    model = ols(f'{result_column} ~ {variable_one} + {variable_2} + {variable_one}:{variable_2}',
                data=dataframe).fit()
    results = sm.stats.anova_lm(model, typ=2)
    print(results)
    # Perform Tukey's HSD labeled_schema_aln.csv for multiple comparisons
    mc_results_one = mc.pairwise_tukeyhsd(dataframe[result_column], dataframe[variable_one])
    print(mc_results_one)
    mc_results_two = mc.pairwise_tukeyhsd(dataframe[result_column], dataframe[variable_2])
    print(mc_results_two)

def nanosight_csv_to_dataframe(summary_csv:str):
    dataframe = pd.read_csv(summary_csv, skiprows=78)
    dataframe = dataframe.drop(
        columns=['Standard Error', 'Unnamed: 8', 'Concentration (particles / ml)', 'Concentration (particles / ml).1',
                 'Concentration (particles / ml).2', 'Concentration (particles / ml).3',
                 'Concentration (particles / ml).4'])
    dataframe = dataframe.dropna()
    dataframe['Bin centre (nm)']=dataframe['Bin centre (nm)'].astype(float)
    aggregation_sum=dataframe[dataframe['Bin centre (nm)'] > nanosight_size_cutoff]['Concentration average'].sum()
    return dataframe,aggregation_sum

def nanosight_aggregation_sum(summary_files:list):
    for file in summary_files:
        virus=file.split()[0]
        virus_data,aggregate_sum=nanosight_csv_to_dataframe(file)
        max_size=virus_data[virus_data['Concentration average']>0]['Bin centre (nm)'].max()
        print(virus,aggregate_sum,max_size)

#Antibody fluorescence analysis
#343.5nm nanosight SARS 90th percentile

if __name__ == '__main__':
    particle_data = pd.read_csv(r"C:\Users\jamel\OneDrive\Documents\Stats Meeting\v200globe_total_bulk_particles.csv")
    fov_intensities = pd.read_csv(r"C:\Research\immunofluorescence\fov_antibody_intensity.csv")
    antibody_intensities = pd.read_csv(r"C:\Research\immunofluorescence\0_5_percent_threshold_antibody_particles.csv")
    nanosight_size_cutoff = 343.5
    # labeled_schema_aln.csv for variance homogeneity and it passed
    # virus_areas=[df for label,df in particle_data.groupby('Virus')['Area']]
    # w,alp=levene(*virus_areas,center='median')
    # for csv in os.listdir('room_temp_nanosight'):

    # nanosight_data.plot(x='Bin centre (nm)',y='Concentration average')
    # plt.xticks(range(0,900,100),range(0,900,100))
    # plt.show()

    # for virus,virus_data in antibody_intensities.groupby('Virus'):
    #     virus_data.plot(kind='hist',y='IntensityBackSub',title=virus)
    #     plt.show()
    filtered_inten=filter_temp_from_particle_table(antibody_intensities,25)
    # print(filtered_inten['IntensityBackSub'].sum())
    print(filtered_inten.count())
    # DLS_Analysis.Read_Histogram.overlapping_line_plot(antibody_intensities,'Temperature','IntensityBackSub','Virus')
    # fov_intensities=fov_intensities[fov_intensities['Temperature'] == 25]
    # print(fov_intensities.groupby("Virus")['BG_Sub_Inten'].mean())
    # groupby("col2").plot(kind="bar", title="DataFrameGroupBy Plot")
    # two_way_anova(particle_data, 'Virus', 'Temperature', 'Area')
    #look at just 25-75
    # virus_tables = filter_temp_from_particle_table(particle_data, '25cel')
    # for virus_data in virus_tables:
    #     print(len(virus_data[virus_data['Area'] <= 31]) / len(virus_data))
