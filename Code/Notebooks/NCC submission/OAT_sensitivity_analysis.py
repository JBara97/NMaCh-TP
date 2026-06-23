import numpy as np
from util.dynamics import temp_response, environmental_update
from scipy.optimize import fsolve, curve_fit

### BASELINE SIMULATIONS
years = np.arange(1850,1850+len(temp))
no_years = len(years)
# Modifying collapse probabilities for each year according to each scenario
p_c = np.ones((no_years,len(df_TEs)))
p_c *= q_c

p_r = np.ones((no_years,len(df_TEs)))
p_r *= q_r

for t in range(no_years):
    p_c[t,:] *= temp_response(temp[t], T_lim, betas)
    p_r[t,:] *= (1-temp_response(temp[t], T_lim, betas))

D_baseline = np.zeros((len_sample,len(df_TEs),no_years))
for r in range(len_sample):
    for t in range(1,no_years):
        D_baseline[r,:,t] = environmental_update(D_baseline[r,:,t-1],A[r],pc=p_c[t,:],pr=p_r[t,:])

### SENSITIVITY ANALYSIS
## Burn in with historical data (up to and inc. 2015)
years_burn_in = np.arange(1850,2016) # up to and including 2015
no_burn_years = len(years_burn_in)
pc_burn = np.ones((no_burn_years,len(df_TEs)))
pc_burn *= q_c
pr_burn = np.ones((no_burn_years,len(df_TEs)))
pr_burn *= q_r
for t in years_burn_in-1850:
    pc_burn[t,:] *= temp_response(temp[t], T_lim, betas)
    pr_burn[t,:] *= (1-temp_response(temp[t], T_lim, betas))

D_burn = np.zeros((len_sample,len(df_TEs)))
D_burn_noint = np.zeros(len(df_TEs))

for r in range(len_sample):
    for t in range(1,no_burn_years):
        D_burn[r,:] = environmental_update(D_burn[r,:],
                                           A[r],
                                           pc=pc_burn[t,:],
                                           pr=pr_burn[t,:])
for t in range(1,no_burn_years):
    D_burn_noint = environmental_update(D_burn_noint,
                                        np.zeros((len(df_TEs),len(df_TEs))),
                                        pc=pc_burn[t,:],
                                        pr=pr_burn[t,:])

## Simulation for 2015 onwards
years = np.arange(2015,1850+len(temp))
no_years=len(years)
D = [np.zeros((len(pars[i]),len_sample,len(df_TEs),no_years)) for i in range(len(popt))]
D_noint = [np.zeros((len(pars[i]),len(df_TEs),no_years)) for i in range(len(popt))]

for i in range(len(popt)):
    for r in range(len_sample):
        D[i][:,r,:,0] = D_burn[r,:]
    D_noint[i][:,:,0] = D_burn_noint[:]

for i in tqdm.tqdm(range(len(popt))):
    for j in range(len(pars[i])):
        for t in range(1,no_years):
            D_noint[i][j,:,t] = environmental_update(D_noint[i][j,:,t-1],
                                                     np.zeros((len(df_TEs),len(df_TEs))),
                                                     pc=pcs[i][j,t,:],pr=prs[i][j,t,:])
            for r in range(len_sample):
                D[i][j,r,:,t] = environmental_update(D[i][j,r,:,t-1],A[r],
                                                     pc=pcs[i][j,t,:],pr=prs[i][j,t,:])