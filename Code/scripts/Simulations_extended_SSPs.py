import numpy as np
import pandas as pd
import io
import scipy as sp

import sys
sys.path.insert(1, '/home/jbara/Code/NMaCh-TP/Code/')
from util.dynamics import environmental_update
from util.dynamics import temp_response

# --------- READING DATA --------------
# Reading data about individual tipping elements:
folder_data = '../../Data/'
df_TEs = pd.read_excel(folder+'TEsummary.xlsx',header=[0],index_col=0)
T_lim = df_TEs['Threshold (°C) – Est.'].to_numpy() # Critical threshold for tipping
betas = df_TEs['Beta (1/°C)'].to_numpy() # Logistic parameter (temp. sensitivity)
q_c = 1/(1+df_TEs['Collapse (years)'].to_numpy()) # collapse leverage
p_r = 1/(1+df_TEs['Recovery (years)'].to_numpy())

# Reading extended SSP temperatures
years = np.arange(1850, 2351) 
no_years = len(years)
df = pd.read_excel('../../Data/extendedSSPs.xlsx')
scenarios = df['scenario'].unique()
dfs = {s: df[df['scenario']==s].set_index('quantile').iloc[:,9:].transpose() for s in scenarios}

# Loading pre-saved interaction matrices
len_sample=500
A = [np.loadtxt('../../Data/Interaction matrices/A_%s.csv'%r) for r in range(len_sample)]

# -------- SIMULATIONS --------------
# Modifying collapse probabilities for each year according to each scenario
p_c = np.ones((len(scenarios),no_years,len(df_AM22)))
p_c *= q_c
for t in range(no_years):
    for i,s in enumerate(scenarios):
        p_c[i,t,:] *= temp_response(dfs[s]['mean'].to_numpy()[t], T_lim, betas)

# Simulating dynamics for TEs with interactions
D = np.zeros((len(scenarios),len_sample,len(df_TEs),no_years))
for r in (range(len_sample)):
    for j in range(len(scenarios)):
        for t in range(1,no_years):
            D[j,r,:,t] = environmental_update(D[j,r,:,t-1],A[r],p_c=p_c[j,t,:],p_r=p_r)

# Simulating dynamics for TEs without interactions
D_noint = np.zeros((len(scenarios),len(df_TEs),no_years))
for t in (range(1,no_years)):
    for j in range(len(scenarios)):
        D_noint[j,:,t] = environmental_update(D_noint[j,:,t-1],np.zeros((len(df_TEs),len(df_TEs))),p_c=p_c[j,t,:],p_r=p_r)