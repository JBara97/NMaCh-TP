import numpy as np
import pandas as pd
import io
import scipy as sp
from scipy.stats import qmc

import sys
sys.path.insert(1, '../')
from util import stationary

# --------- READING DATA --------------
# Reading data about individual tipping elements:
df_TEs = pd.read_excel('../../Data/TEsummary.xlsx',header=[0],index_col=0)
T_lim = df_TEs['Threshold (°C) – Est.'].to_numpy() # Critical threshold for tipping
betas = df_TEs['Beta (1/°C)'].to_numpy() # Logistic parameter (temp. sensitivity)
q_c = 1/(1+df_TEs['Collapse (years)'].to_numpy()) # collapse leverage
p_r = 1/(1+df_TEs['Recovery (years)'].to_numpy())

# Loading pre-saved interaction matrices
A = [np.loadtxt('../../Data/Interaction matrices/A_%s.csv'%r) for r in range(500)]

# -------- NUMERICS ------------------
# Numerically solving the coupled equations for equilibrium distribution
T = np.arange(0,10,step=0.001) # Range of fixed temperatures
D_theory = np.array([[stationary.stat_dist_theory(A[r],q_c,p_r,t,T_lim,betas) 
                     for t in T] for r in (range(500))])
D_theory_noint = np.array([stationary.stat_dist_theory(np.zeros(9),q_c,p_r,t,T_lim,betas)
                           for t in (T)])


### OLD VERSION IGNORE
# D_theory = np.array([[sp.optimize.fsolve(stat_dist_theory_res, 
#                                         np.ones(9), 
#                                         args=(A[r],util.logistic(t, T_lim,betas)*q_c,p_r)) 
#                      for t in T] for r in tqdm.tqdm(range(500))])

# D_theory_noint = np.array([sp.optimize.fsolve(stat_dist_theory_res, 
#                                               np.ones(9), 
#                                               args=(np.zeros(9),util.logistic(t, T_lim,betas)*q_c,p_r)) 
#                            for t in tqdm.tqdm(T)])