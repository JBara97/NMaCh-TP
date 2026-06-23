import numpy as np
from pathlib import Path
import pandas as pd
import sys

### PROVIDE scenarios 
# Robin Lamboll, Joeri Rogelj, & Carl-Friedrich Schleussner. (2025). 
    # Scenario emissions and temperature data for PROVIDE project (v1.3.0) [Data set]. 
    # Zenodo. https://doi.org/10.5281/zenodo.11108741
# https://pure.iiasa.ac.at/id/eprint/18130/1/essoar.10511875.2.pdf for description
# https://zenodo.org/records/11108741 for actual data set

folder = '../../../Data/PROVIDE 2025/'
cols = list(pd.read_csv(folder+'post-2100_inc_redundanttier1_temperature_summary.csv',nrows=1))
df = pd.read_csv(folder+'post-2100_inc_redundanttier1_temperature_summary.csv',
                 usecols = [c for c in cols if (c != 'quantile') and (c!= 'atomic_scen')]).set_index('scenario').T
df.index.names=['year']
scenarios = df.columns.values

### GTPR 2025
df_TEs = pd.read_excel('../../../Data/TEsummary.xlsx',header=[0],index_col=0,sheet_name='GTPR25')
df_TEs[['Threshold (°C) – Min.', 
        'Threshold (°C) – Est.', 
        'Threshold (°C) – Max.']] = df_TEs[['Threshold (°C) – Min.', 
                                            'Threshold (°C) – Est.', 
                                            'Threshold (°C) – Max.']].interpolate(method='linear', axis=1, limit_direction='backward')
T_lim = df_TEs['Threshold (°C) – Est.'].to_numpy() # Critical threshold for tipping
betas = df_TEs['Beta (1/°C)'].to_numpy() # Logistic parameter (temp. sensitivity)
q_c = 1/(1+df_TEs['Collapse (years)'].to_numpy()) # collapse leverage
q_r = 1/(1+df_TEs['Recovery (years)'].to_numpy())

### Pre-generated interaction matrices from Wunderling 2024
len_sample=500
A = [np.loadtxt('../../../Data/Interaction matrices/Wunderling 2024/A_%s.csv'%r) for r in range(len_sample)]