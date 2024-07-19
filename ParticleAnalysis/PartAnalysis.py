import os
import pandas as pd
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd

# parent_directory=r"C:\Research\immunofluorescence"
# for root, dirs, files in os.walk(parent_directory):
#     if not dirs:
#         file=next((file for file in files if file.endswith('tif')))
#         old_path=os.path.join(root,file)
#         _,_,_,virus,temperature,fov=root.split('\\')
#         temperature=temperature[0:2]
#         virus=virus.split('_')[-1]
#         new_path=os.path.join(parent_directory,f'{virus}_{temperature}_{fov}.ome.tif')
#         os.system(f'copy {old_path} {new_path}')

# im trying a two way anova
particle_data=pd.read_csv('v200globe_total_bulk_particles.csv')
particle_data['Virus']=particle_data['Virus'].astype('category')
particle_data['Temperature']=particle_data['Temperature'].astype('category')
print(particle_data.dtypes)
# Performing two-way ANOVA
model = ols('Area ~ C(Virus) + C(Temperature) + C(Virus):C(Temperature)', data=particle_data).fit()
results=sm.stats.anova_lm(model, typ=2)
print(results)