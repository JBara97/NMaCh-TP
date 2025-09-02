import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import networkx as nx
import matplotlib
import scipy as sp
from sklearn.linear_model import LogisticRegression
import util

#### ------------ Useful functions --------------------
def min_or_max_to_num(valstr):
    if type(valstr)!=str:
        return valstr
    elif valstr.isnumeric():
        return float(valstr)
    else:
        return np.nan
        
def str_to_est(eststr,minstr,maxstr):
    if np.isnan(eststr):
        return np.exp(np.nanmean(np.log([minstr,maxstr])))
    else:
        return eststr

def sigmoid_beta(Tmin,Test,Tmax,value_at_min=0.25):
    # print([type(Tmin),type(Test),type(Tmax)])
    if type(Test)!=str:
        if type(Tmin)!= str:
            return np.log(1/value_at_min - 1)/(Test-Tmin)
        elif type(Tmax) != str:
            return np.log(1/(1-value_at_min) - 1)/(Test-Tmax)
        else:
            return np.nan
    else:
        return np.nan

def TE_recovery_time(hyst_str, time_est, time_max):
    if hyst_str == 'Y':
        return float(time_est) # For hysteresis, the recovery is the same as collapse time scale
    elif hyst_str == 'M' or hyst_str =='P':
        return float(time_max) # For partial or maybe, the recovery is the upper bound on time scale
    else:
        return 10000 # ethical time horizon of 10 ky (following Armstrong Mckay, who relaxed it from Lenton et al. 2008)

def weight_dist(resp,resp_str):
    strength_map = {'Strong':10,'Moderate':5,'Weak, Moderate':3,'Weak':2,'Unclear':1.5}
    # The map uses as reference values:
        # Strong = 10 - GIS -> AMOC: Strong (Wunderling 2024) and PF = [1,10] (Kriegler 2009)
        # Moderate = 5 - the second highest absolute PF value in Kriegler 2009 is 5 (ENSO -> WAIS)
        # Weak, Moderate = 3 - WAIS -> AMOC: M/W (Wunderling 2024) and PF = [0.3,3] (Kriegler 2009)
        # Unclear = 1.5 - AMOC -> WAIS: Unclear (Wunderling 2024) and PF = [1,1.5] (Kriegler 2009)
    str_val = strength_map[resp_str]
    if resp == 'Unclear':
        return (-1*str_val,str_val)
    else:
        resp_sign_map = {'Destabilizing':+1,'Stabilizing':-1}
        # The minimum absolute value is 0.1 since that's the smallest value in the PF's of Kriegler 2009 (AMOC -> GIS)
        bounds = resp_sign_map[resp]*np.array([0.1,str_val])
        return (min(bounds),max(bounds))

def logistic(T,T0,beta):
    return 1/(1+np.exp(beta*(T0-T)))

#### ------------ Armstrong McKay 2022 Data ------------------------------------------
fpath1 = '../Data/ArmstrongMcKay2022/science.abn7950_data_s1.xlsx'
df_AM22 = pd.read_excel(fpath1,sheet_name='4) Estimates',header=[5,6],index_col=1,nrows=33, engine='openpyxl')
df_AM22=df_AM22.replace('\n','',regex=True)
df_AM22.columns = [' - '.join(col).strip() for col in df_AM22.columns.values]
df_AM22 = df_AM22.rename(columns={'Category - Unnamed: 0_level_1':'Category',
                        'Proposed \nClimate Tipping Element \n[and Tipping Point]* - Unnamed: 2_level_1':'Tipping Element',
                        'Activation Time - Unnamed: 11_level_1':'Activation Time',
                        'Dynamics - Self-\nperpetuate?':'Dynamics - Self-perpetuate?',
                        'Dynamics - Rate-\ndependent?':'Dynamics - Rate-dependent?',
                        'Comments - Unnamed: 15_level_1':'Comments'})
df_AM22.iloc[1:9,0] = 'Global Core Tipping Elements'
df_AM22.iloc[10:16,0] = 'Regional Impact Tipping Elements'
df_AM22.iloc[16:21,0] = 'Threshold-free nonlinear feedbacks'
df_AM22.iloc[21:30,0] = 'Uncertain Potential Tipping Elements'
df_AM22.iloc[30:,0] = 'Unlikely to be a Tipping Element'

# Converting strings of timescales (which sometimes contains for example ky as 1000y) into numerical value
df_AM22[['Timescale (y) - Est.', 'Timescale (y) - Min.', 'Timescale (y) - Max.']] = df_AM22[['Timescale (y) - Est.', 'Timescale (y) - Min.', 'Timescale (y) - Max.']].replace('k','000',regex=True)
df_AM22[['Timescale (y) - Min.', 'Timescale (y) - Max.']] = df_AM22[['Timescale (y) - Min.', 'Timescale (y) - Max.']].replace('?',np.nan)
df_AM22[['Timescale (y) - Min.', 'Timescale (y) - Max.']] = df_AM22[['Timescale (y) - Min.', 'Timescale (y) - Max.']].replace('/',np.nan,regex=True)
df_AM22['Timescale (y) - Min.'] = df_AM22['Timescale (y) - Min.'].apply(min_or_max_to_num)
df_AM22['Timescale (y) - Max.'] = df_AM22['Timescale (y) - Max.'].apply(min_or_max_to_num)
df_AM22['Timescale (y) - Est.'] = df_AM22['Timescale (y) - Est.'].apply(min_or_max_to_num)

estimates = []
for x in df_AM22[['Timescale (y) - Est.', 'Timescale (y) - Min.', 'Timescale (y) - Max.']].values:
    estimates.append(str_to_est(x[0],x[1],x[2]))
df_AM22['Timescale (y) - Est.'] = estimates

rev_nums = []
for x in df_AM22[['Dynamics - Hysteresis/ bistability?','Timescale (y) - Est.', 'Timescale (y) - Max.']].values:
    rev_nums.append(TE_recovery_time(x[0],x[1],x[2]))
df_AM22['RecoveryTime'] = rev_nums
# Since REEF has no upper bounds we conservatively use 5*estimated timescale for its recovery time
df_AM22.loc['REEF','RecoveryTime'] = 5*df_AM22.loc['REEF']['Timescale (y) - Est.'] 

df_AM22['CollapseTime'] = df_AM22['Timescale (y) - Est.']

df_AM22['T_lim'] = df_AM22['Threshold (°C) - Est.']
# Manually setting ENSO's T_lim as the midpt of its bounds
df_AM22.loc['ENSO','T_lim'] = np.mean([df_AM22.loc['ENSO']['Threshold (°C) - Min.'],
                                       df_AM22.loc['ENSO']['Threshold (°C) - Max.']],dtype=float)

# Inferring sigmoid beta parameters
betas = []
value_at_min=0.1
for x in df_AM22[['Threshold (°C) - Min.','Threshold (°C) - Est.','Threshold (°C) - Max.']].values:
    betas.append(float(sigmoid_beta(x[0],x[1],x[2],value_at_min=value_at_min)))
df_AM22['Sigmoid Beta'] = betas
Tmin = df_AM22.loc['ENSO']['Threshold (°C) - Min.']
Tmax = df_AM22.loc['ENSO']['Threshold (°C) - Max.']
df_AM22.loc['ENSO','Sigmoid Beta'] = np.log(1/value_at_min - 1)/(np.mean([Tmin,Tmax])-Tmin)

# Manually adding recovery times from Nico
df_AM22['RecoveryTime (Nico)'] = {'REEF':5000, 'AWSI' : 20, 'PFTP':500, 
                                  'AMOC':50, 'SAHL':50, 'AMAZ':1000,
                                  'ENSO':100, 'WAIS':2e4, 'GrIS':1e5}

# ------------ Wunderling 2024 Data (TE-TE interactions) -------------
fpath2 = '../Data/Wunderling2024/Fig3 edgelist.ods'
df_W24 = pd.read_excel(fpath2,nrows=19,engine='odf')

# Removing edges in which BOTH Response AND (Response Strength are Unclear OR Response Strength are Very Limited)
df_W24 = df_W24[~((df_W24['Response']=='Unclear')&((df_W24['Response Strength']=='Unclear') | (df_W24['Evidence']=='Very Limited')))]

dist_mins = []
dist_maxs = []
for x in df_W24[['Response','Response Strength']].values:
    bounds = weight_dist(x[0],x[1])
    dist_mins.append(bounds[0])
    dist_maxs.append(bounds[1])

df_W24['Dist. Min.'] = dist_mins
df_W24['Dist. Max.'] = dist_maxs

# --------------- Semantic cohesion ---------------
rename_AM22 = {'GrIS':'GIS', # Greenland Ice Sheet
               'PFTP': 'PERM', # (Boreal) Permafrost collapse
               'AWSI':'ASI', # Arctic (Winter) Sea Ice
               'INSM': 'ISM', # Indian Summer Monsoon
               'SAHL': 'WAM' # Sahel & West African Monsoon -> West African Monsoon
              }
rename_W24 = {'Coral':'REEF' # Coral reefs}
             }
df_AM22=df_AM22.rename(index=rename_AM22)
df_W24=df_W24.replace({'Source':rename_W24,'Target':rename_W24})
node_names=list(set(list(df_W24['Source'].unique())+list(df_W24['Target'].unique())))

df_AM22 = df_AM22.loc[list(set(list(df_W24['Source'].unique())+list(df_W24['Target'].unique())))]
df_AM22[['Tipping Element','CollapseTime','RecoveryTime','T_lim','Sigmoid Beta']]
df_AM22 = df_AM22.sort_index()

TEs = df_AM22.index.to_list()
TE_abb = {te:te for te in TEs}
TE_abb['GIS'] = 'GIS-TE'
TE_abb['WAIS'] = 'WAIS-TE'
TE_abb['ASI'] = 'ASI-TE'
TE_abb['REEF'] = 'REEF-TE'

TE_ordered = ['AMOC', 'WAM', 'REEF',
              'ENSO', 'AMAZ','PERM',
              'GIS', 'WAIS', 'ASI']
TE_ind = {te:i for i,te in enumerate(TEs)}

# -------------- Graph generation ---------------
## This needs to be done alongside simulations since we need to sample each edge weight from U(min,max)

# -------------- Global Warming forecasts (SSPs) -------------------
## https://catalogue.ceda.ac.uk/uuid/98af2184e13e4b91893ab72f301790db/
fpath = '../Data/IPCC AR6 Warming Data/tas_global_'
scenarios = ['Historical',
            'SSP1_1_9',
            'SSP1_2_6',
            'SSP2_4_5',
            'SSP3_7_0',
            'SSP5_8_5']
df_SSPs = {s:pd.read_csv(fpath+s+'.csv') for s in scenarios}
